#!/usr/bin/env python3
"""
LN metastasis prediction (N0 vs Npos) - multimodel + enhanced visualization
Feature modes:
- combined: protein + metabolite
- protein: protein only
- metabolite: metabolite only

Key outputs:
1) SVM-RFE feature selection (tool only, not main model)
2) Multi-model training & evaluation + AUC bootstrap 95% CI + ACC(0.5)
3) ROC (4 panels): cohort2 / cohort3 / Early subset cohort2 / Early subset cohort3
4) Feature importance (SVM-RFE coef_)
5) PCA-style scatter: x=SEM predicted prob(Npos), y=PC1 (by N_status and by stage)
6) SEM ROC data export + SEM-only mode comparison (2x2)

Dependencies:
pandas, numpy, scikit-learn, matplotlib, seaborn, xgboost(optional)
"""

import os
import warnings
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    roc_auc_score, accuracy_score,
    precision_score, recall_score, f1_score,
    confusion_matrix
)
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_curve, auc, roc_auc_score, accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import GradientBoostingClassifier, StackingClassifier
from sklearn.neural_network import MLPClassifier
from matplotlib.patches import Patch

warnings.filterwarnings("ignore")

# optional xgboost
try:
    from xgboost import XGBClassifier
except ImportError:
    XGBClassifier = None
    print("⚠️ 未安装 xgboost，将跳过 XGB 模型。")

# plot defaults
plt.rcParams["font.size"] = 12
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["axes.titleweight"] = "bold"
sns.set_style("whitegrid")

def _safe_div(a, b):
    return float(a) / float(b) if b != 0 else np.nan


def plot_confusion_matrix_DE(y_train, y_val, p_train, p_val, out_png):
    from sklearn.metrics import confusion_matrix
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    label_order = ["N0", "Npos"]
    display = ["Non-LN", "LN"]

    yhat_train = np.where(p_train >= 0.5, "Npos", "N0")
    yhat_val = np.where(p_val >= 0.5, "Npos", "N0")

    cm_train = confusion_matrix(y_train, yhat_train, labels=label_order)
    cm_val = confusion_matrix(y_val, yhat_val, labels=label_order)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # ---- Training ----
    sns.heatmap(
        cm_train, annot=True, fmt="d", cmap="Purples",
        xticklabels=display, yticklabels=display,
        annot_kws={"size": 20},  # 热图数字
        cbar_kws={"label": "Freq"}, ax=axes[0]
    )
    axes[0].set_title("Training Dataset", fontsize=22, fontweight="bold")
    axes[0].set_xlabel("Predicted", fontsize=18)
    axes[0].set_ylabel("Actual", fontsize=18)
    axes[0].tick_params(axis='both', labelsize=18)  # 坐标刻度

    # 颜色条刻度
    axes[0].collections[0].colorbar.ax.tick_params(labelsize=16)
    axes[0].collections[0].colorbar.ax.set_ylabel("Freq", fontsize=18)

    # ---- Validation ----
    sns.heatmap(
        cm_val, annot=True, fmt="d", cmap="Oranges",
        xticklabels=display, yticklabels=display,
        annot_kws={"size": 20},
        cbar_kws={"label": "Freq"}, ax=axes[1]
    )
    axes[1].set_title("Validation Dataset", fontsize=22, fontweight="bold")
    axes[1].set_xlabel("Predicted", fontsize=18)
    axes[1].set_ylabel("Actual", fontsize=18)
    axes[1].tick_params(axis='both', labelsize=18)

    axes[1].collections[0].colorbar.ax.tick_params(labelsize=16)
    axes[1].collections[0].colorbar.ax.set_ylabel("Freq", fontsize=18)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()


def bootstrap_auc(y_true, y_score, n_bootstraps=1000, seed=44):
    """bootstrap 计算 AUC 的 95% CI"""
    rng = np.random.RandomState(seed)
    y_true = np.asarray(y_true, dtype=int)
    y_score = np.asarray(y_score, dtype=float)
    n = len(y_true)
    aucs = []

    for _ in range(n_bootstraps):
        idx = rng.randint(0, n, n)
        if len(np.unique(y_true[idx])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[idx], y_score[idx]))

    if len(aucs) == 0:
        return np.nan, np.nan
    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)


class BreastCancerLNAnalysis:
    def __init__(self, output_dir="./breast_cancer_figures_LN", feature_mode="combined", seed=44):
        """
        feature_mode:
          - combined: 蛋白+代谢物
          - protein : 仅蛋白
          - metabolite: 仅代谢物
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        self.feature_mode = feature_mode
        self.seed = seed

        self.coh1_data = None
        self.coh2_data = None

        self.keep_proteins = None
        self.keep_metas = None

        self.X_train = None
        self.X_test1 = None
        self.X_test2 = None
        self.y_train = None
        self.y_test1 = None
        self.y_test2 = None

        self.X_train_scaled = None
        self.X_test1_scaled = None
        self.X_test2_scaled = None

        self.selected_features = None
        self.feature_importance = None
        self.scaler = None

        self.additional_models = {}
        self.model_performance = None

        # 显示用 cohort 名称（你也可以按需改成 cohort3/3）
        self.cohort_display = {"coh2": "Cohort 2", "coh3": "Cohort 3"}

    # ============= N0/Npos mapping =============
    @staticmethod
    def _map_N_status(n):
        """N==0 -> N0; N>0 -> Npos; else NaN"""
        if pd.isna(n):
            return np.nan
        try:
            n_val = int(float(n))
        except Exception:
            return np.nan
        return "N0" if n_val == 0 else "Npos"

    @staticmethod
    def _is_early_stage(stage):
        """Early stage definition (string tolerant): {0, IA, IB}"""
        if pd.isna(stage):
            return False
        s = str(stage).strip().upper()
        return s in ["0", "IA", "IB"]

    # ==================== 1) load data ====================
    def load_data(self):
        print("正在加载数据...")

        pro_coh1 = pd.read_csv("pro_coh1.csv")
        pro_coh2 = pd.read_csv("pro_coh2.csv")

        meta_coh1 = pd.read_csv("Meta_coh1.csv")
        meta_coh2 = pd.read_csv("Meta_coh2.csv")

        clinical_coh1 = pd.read_csv("Clinical_coh1.csv")
        clinical_coh2 = pd.read_csv("Clinical_coh2.csv")

        # keep lists (LN version)
        with open("combine_protein1.txt", "r") as f:
            self.keep_proteins = [line.strip() for line in f if line.strip()]

        with open("metabolites1.txt", "r") as f:
            self.keep_metas = [line.strip() for line in f if line.strip()]

        protein_id_col = pro_coh1.columns[0]
        meta_id_col = meta_coh1.columns[0]

        # filter
        pro_coh1 = pro_coh1[pro_coh1[protein_id_col].isin(self.keep_proteins)]
        pro_coh2 = pro_coh2[pro_coh2[protein_id_col].isin(self.keep_proteins)]
        meta_coh1 = meta_coh1[meta_coh1[meta_id_col].isin(self.keep_metas)]
        meta_coh2 = meta_coh2[meta_coh2[meta_id_col].isin(self.keep_metas)]

        # transpose: rows=Sample, cols=Feature
        pro_coh1_t = pro_coh1.set_index(protein_id_col).T
        pro_coh2_t = pro_coh2.set_index(protein_id_col).T
        meta_coh1_t = meta_coh1.set_index(meta_id_col).T
        meta_coh2_t = meta_coh2.set_index(meta_id_col).T

        # scale separately (same style as your BC script)
        scaler_pro = StandardScaler()
        pro1_scaled = pd.DataFrame(scaler_pro.fit_transform(pro_coh1_t),
                                   index=pro_coh1_t.index, columns=pro_coh1_t.columns)
        pro2_scaled = pd.DataFrame(scaler_pro.transform(pro_coh2_t),
                                   index=pro_coh2_t.index, columns=pro_coh2_t.columns)

        scaler_meta = StandardScaler()
        meta1_scaled = pd.DataFrame(scaler_meta.fit_transform(meta_coh1_t),
                                    index=meta_coh1_t.index, columns=meta_coh1_t.columns)
        meta2_scaled = pd.DataFrame(scaler_meta.transform(meta_coh2_t),
                                    index=meta_coh2_t.index, columns=meta_coh2_t.columns)

        # join to clinical by Sample
        coh1_tmp = clinical_coh1.set_index("Sample").join(pro1_scaled, how="inner").join(meta1_scaled, how="inner")
        coh2_tmp = clinical_coh2.set_index("Sample").join(pro2_scaled, how="inner").join(meta2_scaled, how="inner")

        self.coh1_data = coh1_tmp.drop_duplicates()
        self.coh2_data = coh2_tmp.drop_duplicates()

        print(f"{self.cohort_display['coh2']}: {self.coh1_data.shape[0]} 样本, {self.coh1_data.shape[1]} 列")
        print(f"{self.cohort_display['coh3']}: {self.coh2_data.shape[0]} 样本, {self.coh2_data.shape[1]} 列")

        return self.coh1_data, self.coh2_data

    # ==================== 2) prepare features (Malignant only; label=N0/Npos) ====================
    def prepare_features(self):
        print("\n准备特征和标签（N0 vs Npos）...")
        print(f"当前特征模式: {self.feature_mode}")

        clinical_cols = [
            "Label", "Clinical stage", "Molecular Classfication",
            "T", "N", "M", "Age", "BMI", "Cohort"
        ]

        # Malignant only
        c1 = self.coh1_data[self.coh1_data["Label"] == "Malignant"].copy()
        c2 = self.coh2_data[self.coh2_data["Label"] == "Malignant"].copy()

        # N_status
        c1["N_status"] = c1["N"].apply(self._map_N_status)
        c2["N_status"] = c2["N"].apply(self._map_N_status)
        c1 = c1[~c1["N_status"].isna()].copy()
        c2 = c2[~c2["N_status"].isna()].copy()

        print(f"{self.cohort_display['coh2']} N_status 分布：")
        print(c1["N_status"].value_counts())
        print(f"{self.cohort_display['coh3']} N_status 分布：")
        print(c2["N_status"].value_counts())

        # feature columns
        all_feature_cols = [col for col in c1.columns if col not in clinical_cols + ["N_status"]]
        protein_cols = [c for c in all_feature_cols if (self.keep_proteins is not None and c in self.keep_proteins)]
        metabolite_cols = [c for c in all_feature_cols if (self.keep_metas is not None and c in self.keep_metas)]

        if self.feature_mode == "protein":
            feature_cols = protein_cols
        elif self.feature_mode == "metabolite":
            feature_cols = metabolite_cols
        else:
            feature_cols = sorted(list(set(protein_cols + metabolite_cols)))

        if len(feature_cols) == 0:
            raise ValueError(f"{self.feature_mode} 模式下没有可用特征，请检查列名与 keep 列表是否一致。")

        print(f"蛋白特征数: {len(protein_cols)}, 代谢物特征数: {len(metabolite_cols)}")
        print(f"当前模式({self.feature_mode})使用特征数: {len(feature_cols)}")

        X1 = c1[feature_cols].fillna(0)
        y1 = c1["N_status"]
        X2 = c2[feature_cols].fillna(0)
        y2 = c2["N_status"]

        # split cohort2 -> train/test1
        X_train, X_test1, y_train, y_test1 = train_test_split(
            X1, y1, test_size=0.3, random_state=self.seed, stratify=y1
        )

        self.coh1_data = c1
        self.coh2_data = c2

        self.X_train, self.X_test1, self.X_test2 = X_train, X_test1, X2
        self.y_train, self.y_test1, self.y_test2 = y_train, y_test1, y2
        self.common_features = feature_cols

        return X_train, X_test1, X2, y_train, y_test1, y2

    # ==================== 3) feature selection (SVM-RFE tool) ====================
    def feature_selection(self, n_features=10):
        print("\n执行特征选择（SVM-RFE）...")

        n_features = min(n_features, self.X_train.shape[1])
        if n_features < 1:
            raise ValueError("有效特征数小于 1，无法进行 SVM-RFE。")
        print(f"SVM-RFE 目标特征数: {n_features}")

        scaler_rfe = StandardScaler()
        X_scaled_rfe = scaler_rfe.fit_transform(self.X_train)

        svm_for_rfe = SVC(kernel="linear", random_state=self.seed)
        rfe = RFE(estimator=svm_for_rfe, n_features_to_select=n_features, step=10)
        rfe.fit(X_scaled_rfe, self.y_train)

        rfe_features = self.X_train.columns[rfe.support_]
        self.selected_features = list(rfe_features)

        print(f"SVM-RFE 选中 {len(rfe_features)} 个特征：")
        print(self.selected_features)

        # unified scaling (selected features only)
        Xtr_sel = self.X_train[rfe_features]
        Xte1_sel = self.X_test1[rfe_features]
        Xte2_sel = self.X_test2[rfe_features]

        self.scaler = StandardScaler()
        self.X_train_scaled = self.scaler.fit_transform(Xtr_sel)
        self.X_test1_scaled = self.scaler.transform(Xte1_sel)
        self.X_test2_scaled = self.scaler.transform(Xte2_sel)

        # importance by coef_
        svm_rfe_fitted = rfe.estimator_
        coef = svm_rfe_fitted.coef_.ravel()

        self.feature_importance = (
            pd.DataFrame({"feature": rfe_features, "coefficient": coef})
            .sort_values("coefficient", key=np.abs, ascending=False)
        )

        feat_path = os.path.join(self.output_dir, "selected_features.csv")
        self.feature_importance.to_csv(feat_path, index=False)
        print(f"特征重要性已保存: {feat_path}")

        return self.selected_features, self.feature_importance

    # ==================== 4) multi-model evaluation ====================
    def evaluate_models(self):
        print("\n=== 训练并评估多种机器学习模型（N0 vs Npos）===")

        Xtr, Xte1, Xte2 = self.X_train_scaled, self.X_test1_scaled, self.X_test2_scaled
        ytr, yte1, yte2 = self.y_train, self.y_test1, self.y_test2

        models = {}

        glm = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=3000, random_state=self.seed,
        )
        models["GLM"] = glm

        lr = LogisticRegression(
            penalty="l2", solver="liblinear",
            max_iter=3000, random_state=self.seed,
        )
        models["LR"] = lr

        svm = SVC(kernel="linear", probability=True, random_state=self.seed)
        models["SVM"] = svm

        gbm = GradientBoostingClassifier(
            learning_rate=0.02, n_estimators=200, max_depth=2,
            subsample=0.8, random_state=self.seed,
        )
        models["GBM"] = gbm

        if XGBClassifier is not None:
            xgb = XGBClassifier(
                n_estimators=550, max_depth=3, learning_rate=0.05,
                subsample=0.8, colsample_bytree=0.8,
                objective="binary:logistic", eval_metric="logloss",
                random_state=self.seed, n_jobs=-1,
            )
            models["XGB"] = xgb
        else:
            print("⚠️ 未安装 xgboost，跳过 XGB 模型。")

        mlp = MLPClassifier(
            learning_rate_init=0.001, hidden_layer_sizes=(64, 32),
            activation="relu", alpha=10, solver="adam",
            max_iter=2000, random_state=self.seed,
        )
        models["MLP"] = mlp

        base_estimators = [("glm", glm), ("lr", lr), ("svm", svm), ("mlp", mlp), ("gbm", gbm)]
        sem = StackingClassifier(
            estimators=base_estimators,
            final_estimator=LogisticRegression(
                penalty="l2", solver="liblinear", max_iter=3000, random_state=self.seed
            ),
            passthrough=False, n_jobs=1,
        )
        models["SEM"] = sem

        # binary truth (Npos=1)
        y1_bin = (yte1 == "Npos").astype(int).values
        y2_bin = (yte2 == "Npos").astype(int).values

        results = []

        for name, clf in models.items():
            print(f"\n>>> 训练/评估模型: {name}")
            clf.fit(Xtr, ytr)
            self.additional_models[name] = clf

            classes = clf.classes_
            if "Npos" not in classes:
                raise ValueError(f"模型 {name} classes 中没有 Npos: {classes}")
            pos_idx = list(classes).index("Npos")

            proba1 = clf.predict_proba(Xte1)[:, pos_idx]
            proba2 = clf.predict_proba(Xte2)[:, pos_idx]

            auc1 = roc_auc_score(y1_bin, proba1)
            auc2 = roc_auc_score(y2_bin, proba2)
            ci1_low, ci1_high = bootstrap_auc(y1_bin, proba1, seed=self.seed)
            ci2_low, ci2_high = bootstrap_auc(y2_bin, proba2, seed=self.seed)

            pred1 = np.where(proba1 >= 0.5, "Npos", "N0")
            pred2 = np.where(proba2 >= 0.5, "Npos", "N0")
            acc1 = accuracy_score(yte1.values, pred1)
            acc2 = accuracy_score(yte2.values, pred2)

            print(f"  {self.cohort_display['coh2']} AUC: {auc1:.3f} (95% CI {ci1_low:.3f}-{ci1_high:.3f}), ACC={acc1:.3f}")
            print(f"  {self.cohort_display['coh3']} AUC: {auc2:.3f} (95% CI {ci2_low:.3f}-{ci2_high:.3f}), ACC={acc2:.3f}")

            results.append({
                "Model": name,
                "cohort2_AUC": auc1,
                "cohort2_CI_low": ci1_low,
                "cohort2_CI_high": ci1_high,
                "cohort2_ACC": acc1,
                "cohort2_n": len(yte1),
                "cohort3_AUC": auc2,
                "cohort3_CI_low": ci2_low,
                "cohort3_CI_high": ci2_high,
                "cohort3_ACC": acc2,
                "cohort3_n": len(yte2),
            })

        results_df = pd.DataFrame(results)
        self.model_performance = results_df

        out_csv = os.path.join(self.output_dir, "LN_N0_vs_Npos_ML_model_comparison.csv")
        results_df.to_csv(out_csv, index=False)
        print(f"\n✅ 多模型对比结果已保存到: {out_csv}")

        return results_df

    def evaluate_sem_only_on_cohort2(self, save_tag="default", threshold=0.5):
        """
        SEM-only：训练在 Train，评估在 Test1（也就是你脚本里的 cohort2/test1 panel）
        输出指标：AUC+bootstrap CI, ACC, Precision, Recall, F1, Specificity
        并保存每个k对应的性能文件（单行csv）
        """
        if self.X_train_scaled is None or self.X_test1_scaled is None:
            raise RuntimeError("请先运行 feature_selection() 生成 X_train_scaled / X_test1_scaled。")

        # ---- build SEM（与 evaluate_models 里保持一致）----
        glm = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=3000, random_state=self.seed,
        )
        lr = LogisticRegression(
            penalty="l2", solver="liblinear",
            max_iter=3000, random_state=self.seed,
        )
        svm = SVC(kernel="linear", probability=True, random_state=self.seed)
        gbm = GradientBoostingClassifier(
            learning_rate=0.02, n_estimators=200, max_depth=2,
            subsample=0.8, random_state=self.seed,
        )
        mlp = MLPClassifier(
            learning_rate_init=0.001, hidden_layer_sizes=(64, 32),
            activation="relu", alpha=10, solver="adam",
            max_iter=2000, random_state=self.seed,
        )

        base_estimators = [("glm", glm), ("lr", lr), ("svm", svm), ("mlp", mlp), ("gbm", gbm)]
        sem = StackingClassifier(
            estimators=base_estimators,
            final_estimator=LogisticRegression(
                penalty="l2", solver="liblinear", max_iter=3000, random_state=self.seed
            ),
            passthrough=False, n_jobs=1,
        )

        # ---- fit ----
        sem.fit(self.X_train_scaled, self.y_train)

        # ---- predict proba on Test1 ----
        classes = list(sem.classes_)
        if "Npos" not in classes:
            raise ValueError(f"SEM classes 中没有 Npos: {classes}")
        pos_idx = classes.index("Npos")

        p = sem.predict_proba(self.X_test1_scaled)[:, pos_idx]
        y_true = (self.y_test1 == "Npos").astype(int).values
        y_pred = (p >= threshold).astype(int)

        # ---- metrics ----
        auc_val = roc_auc_score(y_true, p) if len(np.unique(y_true)) >= 2 else np.nan
        ci_low, ci_high = bootstrap_auc(y_true, p, seed=self.seed) if np.isfinite(auc_val) else (np.nan, np.nan)

        acc = accuracy_score(y_true, y_pred)
        prec = precision_score(y_true, y_pred, zero_division=0)
        rec = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)

        # specificity = TN / (TN + FP)
        cm = confusion_matrix(y_true, y_pred, labels=[0, 1])
        tn, fp, fn, tp = cm.ravel()
        spec = _safe_div(tn, (tn + fp))

        perf = pd.DataFrame([{
            "feature_mode": self.feature_mode,
            "n_features": len(self.selected_features) if self.selected_features is not None else np.nan,
            "threshold": threshold,
            "AUC": auc_val,
            "AUC_CI_low": ci_low,
            "AUC_CI_high": ci_high,
            "Accuracy": acc,
            "Precision": prec,
            "Recall": rec,
            "F1": f1,
            "Specificity": spec,
            "n_test": int(len(y_true)),
            "TP": int(tp), "FP": int(fp), "TN": int(tn), "FN": int(fn),
            "save_tag": save_tag
        }])

        # ---- save per-k ----
        out_csv = os.path.join(self.output_dir, f"SEM_only_perf_cohort2_{self.feature_mode}_{save_tag}.csv")
        perf.to_csv(out_csv, index=False)
        print(f"✅ SEM-only (Test1) 性能已保存: {out_csv}")

        # 把这个 SEM 暂存（可选）
        self.additional_models["SEM_only_sweep"] = sem

        return perf

    def sweep_feature_numbers_sem_only(
            self,
            feature_counts=(5, 10, 15, 20, 30, 40, 50),
            threshold=0.5
    ):
        """
        循环不同 RFE 特征数：
          - 每次：RFE选k个特征 -> 训练SEM -> 在 Test1 验证（你当前逻辑的内部验证）
          - 每个k都单独保存性能文件
          - 最后汇总保存一份总表 + 供画图使用
        """
        if not hasattr(self, "X_train") or self.X_train is None:
            raise RuntimeError("请先运行 load_data() + prepare_features() 以生成训练集。")

        feature_counts = [int(x) for x in feature_counts]
        feature_counts = sorted(list(set([x for x in feature_counts if x >= 1])))

        all_rows = []
        for k in feature_counts:
            k = min(k, self.X_train.shape[1])
            print(f"\n=== SEM-only 扫描：RFE 特征数 = {k} ===")

            # 1) RFE选特征 + 标准化
            self.feature_selection(n_features=k)

            # 2) SEM-only 验证（Test1）
            perf_k = self.evaluate_sem_only_on_cohort2(save_tag=f"n{k}", threshold=threshold)
            all_rows.append(perf_k)

        sweep_df = pd.concat(all_rows, ignore_index=True)
        sweep_df = sweep_df.sort_values("n_features")

        out_csv = os.path.join(self.output_dir, f"SEM_feature_number_sweep_{self.feature_mode}.csv")
        sweep_df.to_csv(out_csv, index=False)
        print(f"\n✅ SEM-only 不同特征数汇总结果已保存: {out_csv}")

        self.sem_feature_sweep_df = sweep_df
        return sweep_df

    def plot_sem_performance_vs_feature_number(self, sweep_df=None, ylim_auc=(0.5, 1.02)):
        """
        画 SEM 的性能随特征数变化：
          - 图1：AUC(含CI) vs 特征数
          - 图2：Acc/F1/Recall/Precision/Specificity vs 特征数
        """
        if sweep_df is None:
            sweep_df = getattr(self, "sem_feature_sweep_df", None)
        if sweep_df is None or len(sweep_df) == 0:
            raise RuntimeError("请先运行 sweep_feature_numbers_sem_only() 得到 sweep_df。")

        sweep_df = sweep_df.sort_values("n_features")

        # ---- 图1：AUC vs 特征数（带CI）----
        fig1, ax1 = plt.subplots(1, 1, figsize=(7, 5))
        ax1.plot(sweep_df["n_features"], sweep_df["AUC"], marker="o", linewidth=2.6)
        ax1.fill_between(
            sweep_df["n_features"].values,
            sweep_df["AUC_CI_low"].values,
            sweep_df["AUC_CI_high"].values,
            alpha=0.15
        )
        ax1.set_xlabel("# Selected features (SVM-RFE)", fontweight="bold", fontsize=14)
        ax1.set_ylabel("AUC", fontweight="bold", fontsize=14)
        ax1.set_title(
            f"SEM performance vs #Features – Test1 ({self.feature_mode})",
            fontweight="bold", fontsize=14
        )
        ax1.set_ylim(*ylim_auc)
        ax1.grid(True, alpha=0.3)
        plt.tight_layout()

        out1 = os.path.join(self.output_dir, f"SEM_AUC_vs_feature_number_{self.feature_mode}.png")
        fig1.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()
        print(f"✅ 已保存: {out1}")

        # ---- 图2：分类指标 vs 特征数 ----
        fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
        for metric in ["Accuracy", "F1", "Recall", "Precision", "Specificity"]:
            if metric in sweep_df.columns:
                ax2.plot(sweep_df["n_features"], sweep_df[metric], marker="o", linewidth=2.3, label=metric)

        ax2.set_xlabel("# Selected features (SVM-RFE)", fontweight="bold", fontsize=14)
        ax2.set_ylabel("Score", fontweight="bold", fontsize=14)
        ax2.set_title(
            f"SEM metrics vs #Features – Test1 ({self.feature_mode})",
            fontweight="bold", fontsize=14
        )
        ax2.set_ylim(0.0, 1.02)
        ax2.legend(loc="lower right")
        ax2.grid(True, alpha=0.3)
        plt.tight_layout()

        out2 = os.path.join(self.output_dir, f"SEM_metrics_vs_feature_number_{self.feature_mode}.png")
        fig2.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()
        print(f"✅ 已保存: {out2}")

        return [fig1, fig2]

    # ==================== SEM ROC data export (for mode comparison) ====================
    def export_sem_roc_data(self):
        """
        导出 SEM 模型在 cohort2 / cohort3 / Early 子集上的 ROC 所需数据
        Early 定义：恶性样本内，Clinical stage ∈ {0, IA, IB}
        """
        if "SEM" not in self.additional_models:
            self.evaluate_models()

        sem = self.additional_models["SEM"]

        Xtr, Xte1, Xte2 = self.X_train_scaled, self.X_test1_scaled, self.X_test2_scaled
        ytr, yte1, yte2 = self.y_train, self.y_test1, self.y_test2

        # cohort2 / cohort3 binary truth (Npos=1)
        y1_bin = (yte1 == "Npos").astype(int).values
        y2_bin = (yte2 == "Npos").astype(int).values

        classes = sem.classes_
        pos_idx = list(classes).index("Npos")

        proba1 = sem.predict_proba(Xte1)[:, pos_idx]
        proba2 = sem.predict_proba(Xte2)[:, pos_idx]

        # ---- Early masks ----
        # cohort2 early built on (train+test1) indices, consistent with your BC script
        idx_all_c1 = list(ytr.index) + list(yte1.index)
        y_all_c1 = pd.concat([ytr, yte1])

        early_mask_c1 = []
        early_labels_c1 = []
        for sample, nlab in zip(idx_all_c1, y_all_c1):
            stage = self.coh1_data.loc[sample, "Clinical stage"] if sample in self.coh1_data.index else np.nan
            if self._is_early_stage(stage):
                early_mask_c1.append(True)
                early_labels_c1.append(1 if nlab == "Npos" else 0)
            else:
                early_mask_c1.append(False)
                early_labels_c1.append(0)

        early_mask_c1 = np.array(early_mask_c1, dtype=bool)
        early_labels_c1 = np.array(early_labels_c1, dtype=int)

        proba_all_c1 = np.concatenate([
            sem.predict_proba(Xtr)[:, pos_idx],
            sem.predict_proba(Xte1)[:, pos_idx],
        ])
        y_early_c1 = early_labels_c1[early_mask_c1]
        p_early_c1 = proba_all_c1[early_mask_c1]

        # cohort3 early uses external only
        idx_all_c2 = list(yte2.index)
        early_mask_c2 = []
        early_labels_c2 = []
        for sample, nlab in zip(idx_all_c2, yte2):
            stage = self.coh2_data.loc[sample, "Clinical stage"] if sample in self.coh2_data.index else np.nan
            if self._is_early_stage(stage):
                early_mask_c2.append(True)
                early_labels_c2.append(1 if nlab == "Npos" else 0)
            else:
                early_mask_c2.append(False)
                early_labels_c2.append(0)

        early_mask_c2 = np.array(early_mask_c2, dtype=bool)
        early_labels_c2 = np.array(early_labels_c2, dtype=int)

        y_early_c2 = early_labels_c2[early_mask_c2]
        p_early_c2 = proba2[early_mask_c2]

        return {
            "mode": self.feature_mode,
            "cohort2": {"y": y1_bin, "p": proba1},
            "cohort3": {"y": y2_bin, "p": proba2},
            "early_cohort2": {"y": y_early_c1, "p": p_early_c1},
            "early_cohort3": {"y": y_early_c2, "p": p_early_c2},
        }

    # ==================== 5) ROC curves (4 panels, BC-style) ====================
    def plot_roc_curves(self):
        print("\n绘制ROC曲线（4 panels）...")

        if not self.additional_models:
            print("尚未运行 evaluate_models()，自动调用...")
            self.evaluate_models()

        preferred_order = ["GLM", "LR", "SVM", "GBM", "XGB", "MLP", "SEM"]
        model_names = [m for m in preferred_order if m in self.additional_models]

        color_list = ["#8B4CB8", "#2E8B57", "#FF6B35", "#2b93db", "#FFD700", "#00B0F0", "#7B68EE", "#FF1493"]
        color_map = {name: color_list[i % len(color_list)] for i, name in enumerate(model_names)}

        Xtr, Xte1, Xte2 = self.X_train_scaled, self.X_test1_scaled, self.X_test2_scaled
        ytr, yte1, yte2 = self.y_train, self.y_test1, self.y_test2

        y1_bin = (yte1 == "Npos").astype(int).values
        y2_bin = (yte2 == "Npos").astype(int).values

        fig, axes = plt.subplots(1, 4, figsize=(24, 6))

        # ---- Panel 1: cohort2 ----
        ax1 = axes[0]
        for name in model_names:
            clf = self.additional_models[name]
            pos_idx = list(clf.classes_).index("Npos")
            p1 = clf.predict_proba(Xte1)[:, pos_idx]
            fpr, tpr, _ = roc_curve(y1_bin, p1)
            this_auc = auc(fpr, tpr)
            ci_low, ci_high = bootstrap_auc(y1_bin, p1, seed=self.seed)
            ax1.plot(fpr, tpr, lw=2.5, color=color_map[name],
                     label=f"{name}: AUC={this_auc:.3f} (95% CI {ci_low:.3f}-{ci_high:.3f})")

        ax1.plot([0, 1], [0, 1], "gray", linestyle="--", lw=1.5, alpha=0.8)
        ax1.set_xlim(0, 1)
        ax1.set_ylim(0, 1.05)
        ax1.set_xlabel("False positive rate", fontweight="bold")
        ax1.set_ylabel("True positive rate", fontweight="bold")
        ax1.set_title(f"LN (N0 vs N+) – {self.cohort_display['coh2']} ({self.feature_mode})",
                      fontweight="bold", fontsize=14)
        ax1.legend(loc="lower right", fontsize=9)
        ax1.grid(True, alpha=0.3)

        # ---- Panel 2: cohort3 ----
        ax2 = axes[1]
        for name in model_names:
            clf = self.additional_models[name]
            pos_idx = list(clf.classes_).index("Npos")
            p2 = clf.predict_proba(Xte2)[:, pos_idx]
            fpr, tpr, _ = roc_curve(y2_bin, p2)
            this_auc = auc(fpr, tpr)
            ci_low, ci_high = bootstrap_auc(y2_bin, p2, seed=self.seed)
            ax2.plot(fpr, tpr, lw=2.5, color=color_map[name],
                     label=f"{name}: AUC={this_auc:.3f} (95% CI {ci_low:.3f}-{ci_high:.3f})")

        ax2.plot([0, 1], [0, 1], "gray", linestyle="--", lw=1.5, alpha=0.8)
        ax2.set_xlim(0, 1)
        ax2.set_ylim(0, 1.05)
        ax2.set_xlabel("False positive rate", fontweight="bold")
        ax2.set_ylabel("True positive rate", fontweight="bold")
        ax2.set_title(f"LN (N0 vs N+) – {self.cohort_display['coh3']} ({self.feature_mode})",
                      fontweight="bold", fontsize=14)
        ax2.legend(loc="lower right", fontsize=9)
        ax2.grid(True, alpha=0.3)

        # ---- Build Early masks ----
        # cohort2 early: use (train+test1)
        idx_all_c1 = list(ytr.index) + list(yte1.index)
        y_all_c1 = pd.concat([ytr, yte1])
        early_mask_c1 = []
        early_labels_c1 = []

        for sample, nlab in zip(idx_all_c1, y_all_c1):
            stage = self.coh1_data.loc[sample, "Clinical stage"] if sample in self.coh1_data.index else np.nan
            if self._is_early_stage(stage):
                early_mask_c1.append(True)
                early_labels_c1.append(1 if nlab == "Npos" else 0)
            else:
                early_mask_c1.append(False)
                early_labels_c1.append(0)

        early_mask_c1 = np.array(early_mask_c1, dtype=bool)
        early_labels_c1 = np.array(early_labels_c1, dtype=int)

        # cohort3 early: external only
        idx_all_c2 = list(yte2.index)
        early_mask_c2 = []
        early_labels_c2 = []
        for sample, nlab in zip(idx_all_c2, yte2):
            stage = self.coh2_data.loc[sample, "Clinical stage"] if sample in self.coh2_data.index else np.nan
            if self._is_early_stage(stage):
                early_mask_c2.append(True)
                early_labels_c2.append(1 if nlab == "Npos" else 0)
            else:
                early_mask_c2.append(False)
                early_labels_c2.append(0)

        early_mask_c2 = np.array(early_mask_c2, dtype=bool)
        early_labels_c2 = np.array(early_labels_c2, dtype=int)

        # ---- Panel 3: cohort2 Early ----
        ax3 = axes[2]
        if early_mask_c1.sum() > 0:
            X_all_c1 = np.vstack([Xtr, Xte1])
            for name in model_names:
                clf = self.additional_models[name]
                pos_idx = list(clf.classes_).index("Npos")
                p_all = np.concatenate([clf.predict_proba(Xtr)[:, pos_idx], clf.predict_proba(Xte1)[:, pos_idx]])
                y_e = early_labels_c1[early_mask_c1]
                p_e = p_all[early_mask_c1]
                if len(np.unique(y_e)) < 2:
                    continue
                fpr, tpr, _ = roc_curve(y_e, p_e)
                this_auc = auc(fpr, tpr)
                ci_low, ci_high = bootstrap_auc(y_e, p_e, seed=self.seed)
                ax3.plot(fpr, tpr, lw=2.5, color=color_map[name],
                         label=f"{name}: AUC={this_auc:.3f} (95% CI {ci_low:.3f}-{ci_high:.3f})")

            ax3.plot([0, 1], [0, 1], "gray", linestyle="--", lw=1.8, alpha=0.8)
            ax3.set_xlim(0, 1)
            ax3.set_ylim(0, 1.05)
            ax3.set_xlabel("False positive rate", fontweight="bold")
            ax3.set_ylabel("True positive rate", fontweight="bold")
            ax3.set_title(f"Early LN (0/IA/IB) – {self.cohort_display['coh2']} ({self.feature_mode})",
                          fontweight="bold", fontsize=14)
            ax3.legend(loc="lower right", fontsize=9)
            ax3.grid(True, alpha=0.3)

        # ---- Panel 4: cohort3 Early ----
        ax4 = axes[3]
        if early_mask_c2.sum() > 0:
            for name in model_names:
                clf = self.additional_models[name]
                pos_idx = list(clf.classes_).index("Npos")
                p_all = clf.predict_proba(Xte2)[:, pos_idx]
                y_e = early_labels_c2[early_mask_c2]
                p_e = p_all[early_mask_c2]
                if len(np.unique(y_e)) < 2:
                    continue
                fpr, tpr, _ = roc_curve(y_e, p_e)
                this_auc = auc(fpr, tpr)
                ci_low, ci_high = bootstrap_auc(y_e, p_e, seed=self.seed)
                ax4.plot(fpr, tpr, lw=2.5, color=color_map[name],
                         label=f"{name}: AUC={this_auc:.3f} (95% CI {ci_low:.3f}-{ci_high:.3f})")

            ax4.plot([0, 1], [0, 1], "gray", linestyle="--", lw=1.8, alpha=0.8)
            ax4.set_xlim(0, 1)
            ax4.set_ylim(0, 1.05)
            ax4.set_xlabel("False positive rate", fontweight="bold")
            ax4.set_ylabel("True positive rate", fontweight="bold")
            ax4.set_title(f"Early LN (0/IA/IB) – {self.cohort_display['coh3']} ({self.feature_mode})",
                          fontweight="bold", fontsize=14)
            ax4.legend(loc="lower right", fontsize=9)
            ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        out_png = os.path.join(self.output_dir, "ROC_LN_N0_vs_Npos_multimodel_with_early.png")
        fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()
        print(f"✅ 多模型 ROC（含 Early）已保存: {out_png}")
        return fig

    def export_predictions_per_sample(self, threshold=0.5, only_models=None):
        """
        导出每个样本的预测结果（每个模型一份记录）：
        - Train (Cohort1 train)
        - Test1 (Cohort1 test)
        - External (Cohort2)

        输出：
        1) predictions_per_sample_all_models.csv  (long format)
        2) predictions_per_sample_SEM_wide.csv    (SEM only, wide format)
        """
        if not self.additional_models:
            raise ValueError("尚未训练模型，请先运行 evaluate_models()。")

        # 需要哪些模型
        model_names = list(self.additional_models.keys())
        if only_models is not None:
            model_names = [m for m in model_names if m in set(only_models)]

        # 三个数据块
        blocks = [
            ("Cohort1", "Train", self.X_train_scaled, self.y_train, self.X_train.index),
            ("Cohort1", "Test1", self.X_test1_scaled, self.y_test1, self.X_test1.index),
            ("Cohort2", "External", self.X_test2_scaled, self.y_test2, self.X_test2.index),
        ]

        # 你想附带导出的临床字段（按需增减）
        # 注意：你的 coh1_data / coh2_data 是 set_index("Sample") 后 join 的，所以 index 是 SampleID
        extra_cols = ["Clinical stage", "Molecular Classification", "T", "N", "M", "Age", "BMI"]

        rows = []

        for cohort_name, split_name, Xblk, yblk, idx in blocks:
            # 取对应的临床表（用于附加字段）
            clin_df = self.coh1_data if cohort_name == "Cohort1" else self.coh2_data

            for model in model_names:
                clf = self.additional_models[model]
                classes = list(clf.classes_)
                if "Npos" not in classes:
                    raise ValueError(f"模型 {model} 的 classes 不包含 Npos: {classes}")
                pos_idx = classes.index("Npos")

                proba = clf.predict_proba(Xblk)[:, pos_idx]
                pred = np.where(proba >= threshold, "Npos", "N0")

                for sid, yt, pp, pl in zip(idx, yblk.values, proba, pred):
                    rec = {
                        "SampleID": sid,
                        "Cohort": cohort_name,
                        "Split": split_name,
                        "TrueLabel": yt,
                        "Model": model,
                        "Proba_Npos": float(pp),
                        "PredLabel": pl,
                        "Threshold": threshold,
                    }

                    # 附加临床字段（如果能取到）
                    if sid in clin_df.index:
                        for c in extra_cols:
                            if c in clin_df.columns:
                                rec[c] = clin_df.loc[sid, c]
                    rows.append(rec)

        pred_long = pd.DataFrame(rows)

        out1 = os.path.join(self.output_dir, "predictions_per_sample_all_models.csv")
        pred_long.to_csv(out1, index=False)
        print(f"✅ 每样本预测结果（全模型，long表）已保存: {out1}")

        # ---------- SEM 额外导出一个“宽表” ----------
        if "SEM" in model_names:
            sem_df = pred_long[pred_long["Model"] == "SEM"].copy()

            sem_wide = sem_df.pivot_table(
                index=["SampleID", "Cohort", "Split", "TrueLabel"],
                values=["Proba_Npos", "PredLabel"],
                aggfunc="first"
            ).reset_index()

            sem_wide = sem_wide.rename(columns={
                "Proba_Npos": "SEM_Proba_Npos",
                "PredLabel": "SEM_PredLabel"
            })

            out2 = os.path.join(self.output_dir, "predictions_per_sample_SEM_wide.csv")
            sem_wide.to_csv(out2, index=False)
            print(f"✅ 每样本预测结果（SEM，wide表）已保存: {out2}")

        return pred_long

    # ==================== 6) feature importance ====================
    def plot_feature_importance(self):
        print("\n绘制特征重要性（基于 SVM-RFE 的 coef_）...")

        if self.feature_importance is None or self.feature_importance.empty:
            print("没有可用的 feature_importance，先运行 feature_selection()。")
            return None

        fig, ax = plt.subplots(1, 1, figsize=(10, 8))
        top_features = self.feature_importance.head(10)
        colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(top_features)))

        bars = ax.bar(
            range(len(top_features)),
            np.abs(top_features["coefficient"]),
            color=colors, alpha=0.8, edgecolor="black", linewidth=0.5
        )

        ax.set_xlabel("Features", fontsize=20, fontweight="bold")
        ax.set_ylabel("Contribution (|coef|)", fontsize=20, fontweight="bold")
        ax.set_title(f"Contribution for LN prediction ({self.feature_mode})",
                     fontsize=22, fontweight="bold", pad=20)

        def wrap_label(label, width=10):
            return textwrap.fill(label, width=width)

        simplified_names = [wrap_label(f.split("_")[0] if "_" in f else f, width=10)
                            for f in top_features["feature"]]

        ax.set_xticks(range(len(top_features)))
        ax.set_xticklabels(simplified_names, rotation=45, ha="right", fontsize=16)
        ax.tick_params(axis="both", which="major", labelsize=18, width=1.5, length=6)

        for bar, value in zip(bars, np.abs(top_features["coefficient"])):
            ax.text(bar.get_x() + bar.get_width() / 2.0,
                    bar.get_height() + 0.001,
                    f"{value:.3f}",
                    ha="center", va="bottom",
                    fontsize=18, fontweight="bold")

        ax.grid(True, alpha=0.3, axis="y")
        ax.set_ylim(0, max(np.abs(top_features["coefficient"])) * 1.15)

        plt.tight_layout()
        out_png = os.path.join(self.output_dir, "feature_importance_SVMRFE.png")
        fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()
        print(f"✅ 特征重要性图已保存: {out_png}")
        return fig

    # ==================== PCA-style visualization (x=SEM prob(Npos), y=PC1) ====================
    def plot_pca_analysis(self):
        """
        x = SEM predicted probability of Npos
        y = PC1 from PCA on selected features
        PCA fit on (Train + Test1) then project to External (cohort3)
        """
        print("\n绘制PCA分析（x=SEM预测Npos概率, y=PC1）...")

        if not self.additional_models:
            self.evaluate_models()
        if "SEM" not in self.additional_models:
            raise ValueError("additional_models 中没有 SEM，请检查 evaluate_models().")

        sem = self.additional_models["SEM"]
        pos_idx = list(sem.classes_).index("Npos")

        p_train = sem.predict_proba(self.X_train_scaled)[:, pos_idx]
        p_test1 = sem.predict_proba(self.X_test1_scaled)[:, pos_idx]
        p_test2 = sem.predict_proba(self.X_test2_scaled)[:, pos_idx]

        # PCA fit on train+test1
        X_all = np.vstack([self.X_train_scaled, self.X_test1_scaled])
        pca = PCA(n_components=2, random_state=self.seed)
        X_pca = pca.fit_transform(X_all)
        X_train_pca = X_pca[:len(self.X_train_scaled)]
        X_test1_pca = X_pca[len(self.X_train_scaled):]
        X_test2_pca = pca.transform(self.X_test2_scaled)

        # colors by N_status
        n0_color = "#2E8B57"
        npos_color = "#8B4CB8"
        legend_elements = [
            Patch(facecolor=n0_color, label="N0"),
            Patch(facecolor=npos_color, label="N+"),
        ]

        # ---- (1) Train+Val by N_status ----
        fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))
        y_train_colors = [npos_color if lab == "Npos" else n0_color for lab in self.y_train]
        ax1.scatter(p_train, X_train_pca[:, 0], c=y_train_colors, alpha=0.75, s=60,
                    edgecolors="black", linewidth=0.5)
        ax1.axvline(x=0.5, color="gray", linestyle="--", alpha=0.7, linewidth=2)
        ax1.set_xlabel("Predicted value for LN (SEM, P[N+])", fontweight="bold", fontsize=18)
        ax1.set_ylabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontweight="bold", fontsize=18)
        ax1.set_title(f"Train+Validation ({self.cohort_display['coh2']})", fontweight="bold", fontsize=20)
        ax1.legend(handles=legend_elements, loc="upper right")
        ax1.grid(True, alpha=0.3)
        ax1.tick_params(axis="both", which="major", labelsize=14, width=1.5, length=6)
        plt.tight_layout()
        path1 = os.path.join(self.output_dir, "pca_train_labels_Nstatus.png")
        fig1.savefig(path1, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        # ---- (2) External by N_status ----
        fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
        y_test2_colors = [npos_color if lab == "Npos" else n0_color for lab in self.y_test2]
        ax2.scatter(p_test2, X_test2_pca[:, 0], c=y_test2_colors, alpha=0.75, s=60,
                    edgecolors="black", linewidth=0.5)
        ax2.axvline(x=0.5, color="gray", linestyle="--", alpha=0.7, linewidth=2)
        ax2.set_xlabel("Predicted value for LN (SEM, P[N+])", fontweight="bold", fontsize=18)
        ax2.set_ylabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", fontweight="bold", fontsize=18)
        ax2.set_title(f"External ({self.cohort_display['coh3']})", fontweight="bold", fontsize=20)
        ax2.legend(handles=legend_elements, loc="upper right")
        ax2.grid(True, alpha=0.3)
        ax2.tick_params(axis="both", which="major", labelsize=14, width=1.5, length=6)
        plt.tight_layout()
        path2 = os.path.join(self.output_dir, "pca_external_labels_Nstatus.png")
        fig2.savefig(path2, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()


        print("✅ PCA + predicted value 图已保存：")

        return [fig1, fig2]

    # ==================== 7) main run ====================
    def run_complete_analysis(self, n_features=6):
        print(f"\n开始完整分析 - LN N0 vs Npos | 模式: {self.feature_mode}")
        print("=" * 70)

        self.load_data()
        self.prepare_features()
        self.feature_selection(n_features=n_features)

        perf_df = self.evaluate_models()
        self.plot_roc_curves()
        self.plot_feature_importance()
        self.plot_pca_analysis()

        # export SEM ROC data (for mode-level summary)
        sem_roc_data = self.export_sem_roc_data()

        print("\n=== 各模型性能列表 ===")
        print(perf_df)

        print("\n=== 分析完成 ===")
        print(f"所有结果已保存到: {self.output_dir}/")

        self.export_predictions_per_sample(threshold=0.5)

        # ===== 自动画 D/E 混淆矩阵（SEM）=====
        sem = self.additional_models["SEM"]
        pos_idx = list(sem.classes_).index("Npos")

        p_train = sem.predict_proba(self.X_train_scaled)[:, pos_idx]
        p_val = sem.predict_proba(self.X_test1_scaled)[:, pos_idx]

        plot_confusion_matrix_DE(
            y_train=self.y_train.values,
            y_val=self.y_test1.values,
            p_train=p_train,
            p_val=p_val,
            out_png=os.path.join(self.output_dir, "ConfusionMatrix_DE.png")
        )

        return {
            "selected_features": self.selected_features,
            "feature_importance": self.feature_importance,
            "model_performance": perf_df,
            "sem_roc_data": sem_roc_data,
        }


def plot_sem_mode_roc_comparison(results_all, outdir="./breast_cancer_figures_LN_SEM_compare"):
    """
    SEM-only, compare modes (combined/protein/metabolite) in 2×2 panels:
    cohort2 / cohort3 / early_cohort2 / early_cohort3
    """
    os.makedirs(outdir, exist_ok=True)

    sem_data = {m: results_all[m]["sem_roc_data"] for m in ["combined", "protein", "metabolite"]}
    colors = {"combined": "#8B4CB8", "protein": "#2b93db", "metabolite": "#FF6B35"}

    panels = [
        ("cohort2", "LN (N0 vs N+) – Cohort 2"),
        ("cohort3", "LN (N0 vs N+) – Cohort 3"),
    ]
    mode_display = {
        "combined": "ProMeta-BC-LN",
        "protein": "Protein-only",
        "metabolite": "Metabolite-only"
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes = axes.flatten()

    for ax, (key, title) in zip(axes, panels):
        for mode in ["combined", "protein", "metabolite"]:
            y = sem_data[mode][key]["y"]
            p = sem_data[mode][key]["p"]

            if len(y) == 0 or len(np.unique(y)) < 2:
                continue

            fpr, tpr, _ = roc_curve(y, p)
            this_auc = auc(fpr, tpr)
            ci_low, ci_high = bootstrap_auc(y, p)

            ax.plot(
                fpr, tpr, lw=2.8, color=colors[mode],
                label=f"{mode_display[mode]}: AUC={this_auc:.3f} \n(95% CI {ci_low:.3f}-{ci_high:.3f})"
            )

        ax.plot([0, 1], [0, 1], "gray", linestyle="--", lw=1.6, alpha=0.8)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("False positive rate", fontweight="bold", fontsize=18)
        ax.set_ylabel("True positive rate", fontweight="bold", fontsize=18)
        ax.tick_params(axis="both", which="major", labelsize=14, width=1.5, length=6)
        ax.set_title(title, fontweight="bold", fontsize=20)
        ax.legend(loc="lower right", fontsize=16)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_path = os.path.join(outdir, "SEM_only_ROC_compare_modes.png")
    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    print(f"✅ SEM-only 三模式 ROC（2×2）对比图已保存: {save_path}")
    return fig




def main():
    results_all = {}
    feature_counts = (2,3,4,5, 6, 7, 8, 9, 10, 11, 12, 13)

    for mode in ["combined", "protein", "metabolite"]:
        outdir = f"./breast_cancer_figures_LN_{mode}"
        analyzer = BreastCancerLNAnalysis(output_dir=outdir, feature_mode=mode, seed=44)

        # 你原来的完整流程（不动）
        results_all[mode] = analyzer.run_complete_analysis(n_features=6)

        # ✅ SEM-only sweep（放 main 里）
        sweep_df = analyzer.sweep_feature_numbers_sem_only(
            feature_counts=feature_counts,
            threshold=0.5
        )
        analyzer.plot_sem_performance_vs_feature_number(sweep_df=sweep_df)

    # 你原来的 SEM 模式对比（不动）
    plot_sem_mode_roc_comparison(results_all, outdir="./breast_cancer_figures_LN_SEM_compare")
    return results_all



if __name__ == "__main__":
    results = main()
