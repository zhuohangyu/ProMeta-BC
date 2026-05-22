#!/usr/bin/env python3
"""
乳腺癌蛋白+代谢物特征预测模型 - 单独保存子图版本
支持三种特征模式：
- combined: 已标准化合并的蛋白 + 代谢物
- protein: 只用蛋白特征
- metabolite: 只用代谢物特征

不再设定“主模型”，仅做多模型性能比较。
"""
import textwrap
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import (
    GradientBoostingClassifier,
    StackingClassifier
)
from sklearn.metrics import (
    accuracy_score, f1_score, recall_score, precision_score,
    confusion_matrix, roc_auc_score
)
from sklearn.neural_network import MLPClassifier
from matplotlib.patches import Patch
import os
import warnings

warnings.filterwarnings('ignore')

# 尝试导入 XGBoost
try:
    from xgboost import XGBClassifier
except ImportError:
    XGBClassifier = None
    print("⚠️ 未安装 xgboost，将跳过 XGB 模型。")

# 设置绘图参数
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
sns.set_style("whitegrid")


def bootstrap_auc(y_true, y_pred, n_bootstraps=1000, seed=40):
    """通过 bootstrap 计算 AUC 95% CI"""
    rng = np.random.RandomState(seed)
    aucs = []
    n = len(y_true)

    for i in range(n_bootstraps):
        indices = rng.randint(0, n, n)
        if len(np.unique(y_true[indices])) < 2:
            continue
        aucs.append(roc_auc_score(y_true[indices], y_pred[indices]))

    return np.percentile(aucs, 2.5), np.percentile(aucs, 97.5)



class BreastCancerAnalysis:
    def __init__(self, output_dir="./figures", feature_mode="combined"):
        """
        feature_mode:
          - "combined": 蛋白 + 代谢物
          - "protein": 仅蛋白
          - "metabolite": 仅代谢物
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)

        self.feature_mode = feature_mode   # 特征模式

        # 多模型容器
        self.additional_models = {}
        self.model_performance = None

        # 在 load_data 中会填充
        self.keep_proteins = None
        self.keep_metas = None
        # 显示用 cohort 名称：把原 coh1/coh2 映射成 cohort2/cohort3
        self.cohort_display = {"coh1": "Cohort 2", "coh2": "Cohort 3", "coh3": "Cohort 4 (n=9)"}

    # ==================== 1. 数据加载 ====================
    def load_data(self):
        """加载和预处理数据（蛋白 + 代谢物 + 临床）"""
        print("正在加载数据...")

        # 读取数据文件
        pro_coh1 = pd.read_csv('pro_coh1.csv')
        pro_coh2 = pd.read_csv('pro_coh2.csv')
        pro_coh3 = pd.read_csv('pro_coh3.csv')

        meta_coh1 = pd.read_csv('Meta_coh1.csv')
        meta_coh2 = pd.read_csv('Meta_coh2.csv')
        meta_coh3 = pd.read_csv('Meta_coh3.csv')

        clinical_coh1 = pd.read_csv('Clinical_coh1.csv')
        clinical_coh2 = pd.read_csv('Clinical_coh2.csv')
        clinical_coh3 = pd.read_csv('Clinical_coh3.csv')

        # 读取要保留的蛋白列表
        with open('combine_protein.txt', 'r') as f:
            keep_proteins = [line.strip() for line in f if line.strip()]
        self.keep_proteins = keep_proteins

        protein_col = pro_coh1.columns[0]

        # 筛选蛋白
        pro_coh1 = pro_coh1[pro_coh1[protein_col].isin(keep_proteins)]
        pro_coh2 = pro_coh2[pro_coh2[protein_col].isin(keep_proteins)]
        pro_coh3 = pro_coh3[pro_coh3[protein_col].isin(keep_proteins)]

        pro_coh1.to_csv('pro_coh1_filtered.csv', index=False)
        pro_coh2.to_csv('pro_coh2_filtered.csv', index=False)

        # 读取要保留的代谢物列表
        with open('cluster_1_metabolites_for_pathway_benign.txt', 'r') as f:
            keep_metas = [line.strip() for line in f if line.strip()]
        self.keep_metas = keep_metas

        meta_col = meta_coh1.columns[0]
        meta_coh1 = meta_coh1[meta_coh1[meta_col].isin(keep_metas)]
        meta_coh2 = meta_coh2[meta_coh2[meta_col].isin(keep_metas)]
        meta_coh3 = meta_coh3[meta_coh3[meta_col].isin(keep_metas)]

        meta_coh1.to_csv('Meta_coh1_filtered.csv', index=False)
        meta_coh2.to_csv('Meta_coh2_filtered.csv', index=False)

        # 转置蛋白（行: Sample, 列: Gene）
        pro_coh1_t = pro_coh1.set_index('Gene').T
        pro_coh2_t = pro_coh2.set_index('Gene').T
        pro_coh3_t = pro_coh3.set_index('Gene').T

        # 转置代谢物（行: Sample, 列: Metabolite）
        meta_coh1_t = meta_coh1.set_index('Metabolite').T
        meta_coh2_t = meta_coh2.set_index('Metabolite').T
        meta_coh3_t = meta_coh3.set_index('Metabolite').T

        # 标准化蛋白
        scaler_protein = StandardScaler()
        pro_coh1_scaled = pd.DataFrame(scaler_protein.fit_transform(pro_coh1_t),
                                       index=pro_coh1_t.index, columns=pro_coh1_t.columns)
        pro_coh2_scaled = pd.DataFrame(scaler_protein.transform(pro_coh2_t),
                                       index=pro_coh2_t.index, columns=pro_coh2_t.columns)
        pro_coh3_scaled = pd.DataFrame(scaler_protein.transform(pro_coh3_t),
                                       index=pro_coh3_t.index, columns=pro_coh3_t.columns)

        # 标准化代谢物
        scaler_metabolite = StandardScaler()
        meta_coh1_scaled = pd.DataFrame(scaler_metabolite.fit_transform(meta_coh1_t),
                                        index=meta_coh1_t.index, columns=meta_coh1_t.columns)
        meta_coh2_scaled = pd.DataFrame(scaler_metabolite.transform(meta_coh2_t),
                                        index=meta_coh2_t.index, columns=meta_coh2_t.columns)
        meta_coh3_scaled = pd.DataFrame(scaler_metabolite.transform(meta_coh3_t),
                                        index=meta_coh3_t.index, columns=meta_coh3_t.columns)

        # 合并到临床
        coh1_temp = clinical_coh1.set_index('Sample').join(pro_coh1_scaled, how='inner')
        coh2_temp = clinical_coh2.set_index('Sample').join(pro_coh2_scaled, how='inner')
        coh3_temp = clinical_coh3.set_index('Sample').join(pro_coh3_scaled, how='inner')

        self.coh1_data = coh1_temp.join(meta_coh1_scaled, how='inner')
        self.coh2_data = coh2_temp.join(meta_coh2_scaled, how='inner')
        self.coh3_data = coh3_temp.join(meta_coh3_scaled, how='inner')

        self.protein_scaler = scaler_protein
        self.metabolite_scaler = scaler_metabolite

        print(f"Cohort 1: {self.coh1_data.shape[0]} 样本, {self.coh1_data.shape[1]} 列")
        print(f"Cohort 2: {self.coh2_data.shape[0]} 样本, {self.coh2_data.shape[1]} 列")
        print(f"Cohort 3: {self.coh3_data.shape[0]} 样本, {self.coh3_data.shape[1]} 列")
        print(f"Cohort 1 标签分布:\n{self.coh1_data['Label'].value_counts()}")
        print(f"Cohort 2 标签分布:\n{self.coh2_data['Label'].value_counts()}")
        print(f"Cohort 3 标签分布:\n{self.coh3_data['Label'].value_counts()}")

        return self.coh1_data, self.coh2_data

    # ==================== 2. 特征准备 ====================
    def prepare_features(self):
        """准备特征和标签（根据 feature_mode 切换）"""
        print("\n准备特征和标签...")
        print(f"当前特征模式: {self.feature_mode}")

        clinical_cols = [
            'Label', 'Clinical stage', 'Molecular Classification',
            'T', 'N', 'M', 'Age', 'BMI', 'Cohort'
        ]

        all_feature_cols = [col for col in self.coh1_data.columns if col not in clinical_cols]

        protein_cols = [c for c in all_feature_cols
                        if self.keep_proteins is not None and c in self.keep_proteins]
        metabolite_cols = [c for c in all_feature_cols
                           if self.keep_metas is not None and c in self.keep_metas]

        if self.feature_mode == "protein":
            feature_cols = protein_cols
        elif self.feature_mode == "metabolite":
            feature_cols = metabolite_cols
        else:
            feature_cols = sorted(list(set(protein_cols + metabolite_cols)))

        if len(feature_cols) == 0:
            raise ValueError(
                f"当前特征模式 {self.feature_mode} 下没有可用特征，请检查列名与 keep_proteins/keep_metas 是否一致。"
            )

        print(f"蛋白特征数: {len(protein_cols)}, 代谢物特征数: {len(metabolite_cols)}")
        print(f"当前模式({self.feature_mode})使用特征数: {len(feature_cols)}")

        # coh1: 用来 split train/test1
        X_coh1 = self.coh1_data[feature_cols].fillna(0)
        y_coh1 = self.coh1_data['Label']

        # coh2: 外部验证（原来的 X_test2）
        X_coh2 = self.coh2_data[feature_cols].fillna(0)
        y_coh2 = self.coh2_data['Label']

        # ✅ coh3: 小样本验证（n=9），不算AUC
        X_coh3 = self.coh3_data[feature_cols].fillna(0)
        y_coh3 = self.coh3_data['Label']

        X_train, X_test1, y_train, y_test1 = train_test_split(
            X_coh1, y_coh1, test_size=0.3, random_state=46, stratify=y_coh1
        )

        self.X_train, self.X_test1, self.X_test2, self.X_test3 = X_train, X_test1, X_coh2, X_coh3
        self.y_train, self.y_test1, self.y_test2, self.y_test3 = y_train, y_test1, y_coh2, y_coh3
        self.common_features = feature_cols

        return X_train, X_test1, X_coh2, X_coh3, y_train, y_test1, y_coh2, y_coh3

    # ==================== 3. 特征选择（只作为工具） ====================
    def feature_selection(self, n_features=20):
        """
        使用 SVM-RFE 做特征选择（不把该 SVM 当“主模型”）。
        选出 n_features 个特征，并对选出的特征做标准化。
        """
        print(f"\n执行特征选择（SVM-RFE）...")

        n_features = min(n_features, self.X_train.shape[1])
        if n_features < 1:
            raise ValueError("有效特征数小于 1，无法进行 SVM-RFE。")

        print(f"SVM-RFE 目标特征数: {n_features}")

        # 1) 先标准化
        scaler_rfe = StandardScaler()
        X_scaled_rfe = scaler_rfe.fit_transform(self.X_train)

        # 2) 定义用于 RFE 的基础 SVM
        svm_for_rfe = SVC(kernel='linear', random_state=40)

        # 3) RFE 拟合
        rfe = RFE(estimator=svm_for_rfe, n_features_to_select=n_features, step=10)
        rfe.fit(X_scaled_rfe, self.y_train)

        # 4) 选中的特征
        rfe_features = self.X_train.columns[rfe.support_]
        print(f"SVM-RFE 选中 {len(rfe_features)} 个特征：")
        print(list(rfe_features))

        # 5) 用选中后的特征构建矩阵，并重新标准化，供后续所有模型使用
        X_train_sel = self.X_train[rfe_features]
        X_test1_sel = self.X_test1[rfe_features]
        X_test2_sel = self.X_test2[rfe_features]
        X_test3_sel = self.X_test3[rfe_features]  # ✅新增

        scaler = StandardScaler()
        self.X_train_scaled = scaler.fit_transform(X_train_sel)
        self.X_test1_scaled = scaler.transform(X_test1_sel)
        self.X_test2_scaled = scaler.transform(X_test2_sel)
        self.X_test3_scaled = scaler.transform(X_test3_sel)  # ✅新增
        self.selected_features = list(rfe_features)
        self.scaler = scaler

        # 6) 用 RFE 训练好的 SVM（rfe.estimator_）来拿 coef_
        svm_rfe_fitted = rfe.estimator_
        coef = svm_rfe_fitted.coef_.ravel()

        feature_importance = pd.DataFrame({
            'feature': rfe_features,
            'coefficient': coef
        }).sort_values('coefficient', key=np.abs, ascending=False)
        self.feature_importance = feature_importance

        print("根据 SVM-RFE 估计的特征重要性：")
        print(feature_importance)

        return self.selected_features, self.feature_importance

    def select_model_by_cv_on_training(self, n_splits=5, n_features=20, seed=40):
        """
        只在 training set 上做 5-fold CV 进行“方法选择”。
        每一折都：StandardScaler -> SVM-RFE -> (再StandardScaler) -> fit model -> AUC on val
        返回：cv_results_df（每个模型 mean/std AUC），并设置 self.selected_model_by_cv
        """
        X = self.X_train.copy()
        y = (self.y_train == "Malignant").values.astype(int)

        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)

        # 备选模型（和 evaluate_models 尽量保持一致）
        models = {}
        glm = LogisticRegression(penalty="elasticnet", l1_ratio=0.5, solver="saga",
                                 max_iter=3000, random_state=40)
        lr = LogisticRegression(penalty="l2", solver="liblinear",
                                max_iter=3000, random_state=40)
        svm = SVC(kernel="linear", probability=True, random_state=40)
        gbm = GradientBoostingClassifier(learning_rate=0.01, n_estimators=200,
                                         max_depth=1, subsample=0.8, random_state=40)
        mlp = MLPClassifier(learning_rate_init=0.001, hidden_layer_sizes=(64, 32),
                            activation="relu", alpha=50, solver="adam",
                            max_iter=2000, random_state=40)

        models["GLM"] = glm
        models["LR"] = lr
        models["SVM"] = svm
        models["GBM"] = gbm
        models["MLP"] = mlp

        if XGBClassifier is not None:
            xgb = XGBClassifier(
                n_estimators=300, max_depth=3,
                learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                objective="binary:logistic", eval_metric="logloss",
                random_state=40, n_jobs=-1,
            )
            models["XGB"] = xgb

        base_estimators = [("glm", glm), ("lr", lr), ("svm", svm), ("mlp", mlp), ("gbm", gbm)]
        sem = StackingClassifier(
            estimators=base_estimators,
            final_estimator=LogisticRegression(penalty="l2", solver="liblinear", max_iter=3000, random_state=40),
            passthrough=False,
            n_jobs=1,
        )
        models["SEM"] = sem

        # 存 AUC：dict[model] -> list of fold aucs
        auc_dict = {m: [] for m in models.keys()}
        fold_rows = []  # 新增：存每一折

        for fold, (tr_idx, va_idx) in enumerate(skf.split(X, y), start=1):
            X_tr, X_va = X.iloc[tr_idx], X.iloc[va_idx]
            y_tr, y_va = y[tr_idx], y[va_idx]

            # === (1) fold 内：StandardScaler + SVM-RFE 选特征 ===
            n_sel = min(n_features, X_tr.shape[1])
            scaler_rfe = StandardScaler()
            X_tr_scaled_for_rfe = scaler_rfe.fit_transform(X_tr)

            rfe = RFE(estimator=SVC(kernel="linear", random_state=40),
                      n_features_to_select=n_sel, step=10)
            rfe.fit(X_tr_scaled_for_rfe, y_tr)

            feats = X_tr.columns[rfe.support_]

            # === (2) fold 内：用选中面板再做一次 scaler（训练折 fit）===
            Xtr_sel = X_tr[feats]
            Xva_sel = X_va[feats]
            scaler = StandardScaler()
            Xtr2 = scaler.fit_transform(Xtr_sel)
            Xva2 = scaler.transform(Xva_sel)

            # === (3) fold 内：训练 + 验证 AUC ===
            for name, clf in models.items():
                clf.fit(Xtr2, y_tr)
                proba = clf.predict_proba(Xva2)[:, 1]
                auc_val = roc_auc_score(y_va, proba)
                auc_dict[name].append(float(auc_val))
                fold_rows.append({
                    "Model": name,
                    "Fold": fold,
                    "AUC": float(auc_val),
                    "Feature_Mode": self.feature_mode
                })

        # 汇总
        rows = []
        for name, aucs in auc_dict.items():
            rows.append({
                "Model": name,
                "CV_mean_AUC": float(np.mean(aucs)),
                "CV_std_AUC": float(np.std(aucs, ddof=1)) if len(aucs) > 1 else 0.0,
                "n_folds": len(aucs),
            })
        cv_df = pd.DataFrame(rows).sort_values("CV_mean_AUC", ascending=False)

        # 选最优（也可改成固定 SEM）
        self.selected_model_by_cv = cv_df.iloc[0]["Model"]

        out = os.path.join(self.output_dir, f"CV_model_selection_{self.feature_mode}.csv")
        cv_df.to_csv(out, index=False)
        print(f"✅ 训练集 {n_splits}-fold CV 方法选择结果已保存: {out}")
        print(f"🏆 CV-best model: {self.selected_model_by_cv} | mean AUC={cv_df.iloc[0]['CV_mean_AUC']:.3f}")
        out_csv = os.path.join(
            self.output_dir,
            f"CV5_training_model_selection_{self.feature_mode}.csv"
        )
        cv_df.to_csv(out_csv, index=False)

        print(f"✅ 5-fold CV (training set) model selection results saved to: {out_csv}")

        fold_df = pd.DataFrame(fold_rows)

        out_fold_csv = os.path.join(
            self.output_dir,
            f"CV5_training_foldwise_AUC_{self.feature_mode}.csv"
        )
        fold_df.to_csv(out_fold_csv, index=False)

        print(f"✅ Fold-wise CV AUC saved to: {out_fold_csv}")
        return cv_df

    # ==================== 4. 多模型评估（核心输出） ====================
    def evaluate_models(self, small_cohort_threshold=0.5):
        """
        Cohort1(Test1): AUC + CI
        Cohort2(External): AUC + CI
        Cohort3(Small n=9): 不算AUC，只算正确预测个数（以及Accuracy可选）
        """
        age_rows = []  # 单独存年龄子组结果，不混进 results
        age_subgroup_for = "best"  # "best" 或 "all"
        print("\n=== 训练并评估多种机器学习模型 ===")

        Xtr, Xte1, Xte2 = self.X_train_scaled, self.X_test1_scaled, self.X_test2_scaled
        ytr, yte1, yte2 = self.y_train, self.y_test1, self.y_test2

        # ✅ 小队列
        has_small = hasattr(self, "X_test3_scaled") and self.X_test3_scaled is not None
        if has_small:
            Xte3, yte3 = self.X_test3_scaled, self.y_test3
            y3_bin = (yte3 == "Malignant").values.astype(int)
            n3 = len(y3_bin)
        else:
            Xte3 = yte3 = None
            y3_bin = None
            n3 = 0

        models = {}

        glm = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=3000, random_state=40,
        )
        models["GLM"] = glm

        lr = LogisticRegression(
            penalty="l2", solver="liblinear",
            max_iter=3000, random_state=40,
        )
        models["LR"] = lr

        svm = SVC(kernel='linear', probability=True, random_state=40)
        models["SVM"] = svm

        gbm = GradientBoostingClassifier(
            learning_rate=0.01, n_estimators=200,
            max_depth=1, subsample=0.8, random_state=40
        )
        models["GBM"] = gbm

        if XGBClassifier is not None:
            xgb = XGBClassifier(
                n_estimators=300, max_depth=3,
                learning_rate=0.05, subsample=0.8, colsample_bytree=0.8,
                objective="binary:logistic", eval_metric="logloss",
                random_state=40, n_jobs=-1,
            )
            models["XGB"] = xgb
        else:
            print("⚠️ 未安装 xgboost，跳过 XGB 模型。")

        mlp = MLPClassifier(
            learning_rate_init=0.001,
            hidden_layer_sizes=(64, 32),
            activation="relu",
            alpha=50,
            solver="adam",
            max_iter=2000,
            random_state=40,
        )
        models["MLP"] = mlp

        base_estimators = [
            ("glm", glm), ("lr", lr), ("svm", svm), ("mlp", mlp), ("gbm", gbm),
        ]
        sem = StackingClassifier(
            estimators=base_estimators,
            final_estimator=LogisticRegression(
                penalty="l2", solver="liblinear", max_iter=3000, random_state=40,
            ),
            passthrough=False,
            n_jobs=1,
        )
        models["SEM"] = sem

        results = []

        for name, clf in models.items():
            print(f"\n>>> 训练/评估模型: {name}")
            clf.fit(Xtr, ytr)
            self.additional_models[name] = clf

            # Cohort1(Test1) AUC
            proba1 = clf.predict_proba(Xte1)[:, 1]
            y1_bin = (yte1 == "Malignant").values
            auc1 = roc_auc_score(y1_bin, proba1)
            ci1_low, ci1_high = bootstrap_auc(y1_bin, proba1)

            # Cohort2(External) AUC
            proba2 = clf.predict_proba(Xte2)[:, 1]
            y2_bin = (yte2 == "Malignant").values
            auc2 = roc_auc_score(y2_bin, proba2)
            ci2_low, ci2_high = bootstrap_auc(y2_bin, proba2)

            print(f"  Cohort1 AUC: {auc1:.3f} (95% CI {ci1_low:.3f}-{ci1_high:.3f})")
            print(f"  Cohort2 AUC: {auc2:.3f} (95% CI {ci2_low:.3f}-{ci2_high:.3f})")

            row = {
                "Model": name,
                "Cohort1_AUC": auc1,
                "Cohort1_CI_low": ci1_low,
                "Cohort1_CI_high": ci1_high,
                "Cohort2_AUC": auc2,
                "Cohort2_CI_low": ci2_low,
                "Cohort2_CI_high": ci2_high,
            }

            # ✅ Cohort3 小样本：正确预测个数（不算AUC）
            if has_small and n3 > 0:
                proba3 = clf.predict_proba(Xte3)[:, 1]
                pred3 = (proba3 >= small_cohort_threshold).astype(int)
                correct3 = int((pred3 == y3_bin).sum())
                acc3 = correct3 / n3

                print(f"  Cohort3 Small (n={n3}) correct: {correct3}/{n3} (thr={small_cohort_threshold})")

                row.update({
                    "Cohort3_n": n3,
                    "Cohort3_threshold": small_cohort_threshold,
                    "Cohort3_correct": correct3,
                    "Cohort3_accuracy": acc3,
                })
            else:
                row.update({
                    "Cohort3_n": n3,
                    "Cohort3_threshold": small_cohort_threshold,
                    "Cohort3_correct": np.nan,
                    "Cohort3_accuracy": np.nan,
                })


            results.append(row)

            # --- 导出每个模型在 Cohort1(Test1) / Cohort2(External) 的逐样本预测 ---
            self._export_predictions(
                model_name=name,
                cohort_tag="Cohort1_Test1",
                sample_index=self.y_test1.index,
                y_series=self.y_test1,
                proba=proba1,
                threshold=0.5
            )
            self._export_predictions(
                model_name=name,
                cohort_tag="Cohort2_External",
                sample_index=self.y_test2.index,
                y_series=self.y_test2,
                proba=proba2,
                threshold=0.5
            )

        # ===== 年龄子组统计（coh1_test1 / coh2_external）=====
        do_age = (age_subgroup_for == "all")  # all 模式：每个模型都算
        if do_age:
            age1 = self._get_age_series_for_samples("coh1", self.y_test1.index)
            age2 = self._get_age_series_for_samples("coh2", self.y_test2.index)

            y1 = (self.y_test1 == "Malignant").values.astype(int)
            y2 = (self.y_test2 == "Malignant").values.astype(int)

            for grp_name, mask1, mask2 in [
                ("Age_lt40", (age1 < 40), (age2 < 40)),
                ("Age_ge55", (age1 >= 55), (age2 >= 55)),
            ]:
                # Cohort1 子组
                idx1 = np.where(mask1.values)[0]
                if len(idx1) > 0:
                    m1 = self._calc_binary_metrics(y1[idx1], proba1[idx1], threshold=0.5)
                    results.append({
                        "Model": name,
                        "Subset": grp_name,
                        "Cohort": self.cohort_display["coh1"],
                        "n": int(len(idx1)),
                        "n_pos": int(y1[idx1].sum()),
                        **m1
                    })

                # Cohort2 子组
                idx2 = np.where(mask2.values)[0]
                if len(idx2) > 0:
                    m2 = self._calc_binary_metrics(y2[idx2], proba2[idx2], threshold=0.5)
                    results.append({
                        "Model": name,
                        "Subset": grp_name,
                        "Cohort": self.cohort_display["coh2"],
                        "n": int(len(idx2)),
                        "n_pos": int(y2[idx2].sum()),
                        **m2
                    })

        results_df = pd.DataFrame(results)
        self.model_performance = results_df

        # ========== 选最优模型（默认：Cohort2_AUC 最大） ==========
        # ========== 强制选择 SEM 作为最优模型 ==========
        best_model = "SEM"
        if best_model not in self.additional_models:
            raise RuntimeError("❌ SEM 未训练成功：请确认 evaluate_models() 中 models 包含 SEM 且已 fit。")

        self.best_model_name = best_model

        # 为了打印也能显示对应AUC（如果结果表里有）
        best_metric = "Cohort2_AUC"
        if "Model" in results_df.columns and best_model in set(results_df["Model"]):
            best_row = results_df.loc[results_df["Model"] == best_model].iloc[0]
            val = best_row.get(best_metric, np.nan)
            print(f"\n🏆 Best model (forced): {best_model} | {best_metric}={val:.3f}")
        else:
            print(f"\n🏆 Best model (forced): {best_model}")

        # ========== 再额外导出 best model（固定文件名，便于你后续直接用） ==========
        best_clf = self.additional_models[best_model]

        best_proba1 = best_clf.predict_proba(self.X_test1_scaled)[:, 1]
        best_proba2 = best_clf.predict_proba(self.X_test2_scaled)[:, 1]

        out1 = self._export_predictions(
            model_name=f"BEST_{best_model}",
            cohort_tag="Cohort1_Test1",
            sample_index=self.y_test1.index,
            y_series=self.y_test1,
            proba=best_proba1,
            threshold=0.5
        )
        out2 = self._export_predictions(
            model_name=f"BEST_{best_model}",
            cohort_tag="Cohort2_External",
            sample_index=self.y_test2.index,
            y_series=self.y_test2,
            proba=best_proba2,
            threshold=0.5
        )

        # 如果你也想对 Cohort3 small 导出（可选）
        if hasattr(self, "X_test3_scaled") and self.X_test3_scaled is not None:
            best_proba3 = best_clf.predict_proba(self.X_test3_scaled)[:, 1]
            self._export_predictions(
                model_name=f"BEST_{best_model}",
                cohort_tag="Cohort3_Small",
                sample_index=self.y_test3.index,
                y_series=self.y_test3,
                proba=best_proba3,
                threshold=small_cohort_threshold
            )


        cmp_path = os.path.join(self.output_dir, "ML_model_comparison.csv")
        results_df.to_csv(cmp_path, index=False)
        # ====== 只对 SEM 做年龄子组统计（唯一入口） ======
        best_model = "SEM"
        self.best_model_name = best_model
        best_clf = self.additional_models[best_model]

        # 计算 best model 的 proba（与导出一致）
        best_proba1 = best_clf.predict_proba(self.X_test1_scaled)[:, 1]
        best_proba2 = best_clf.predict_proba(self.X_test2_scaled)[:, 1]

        # 年龄
        age1 = self._get_age_series_for_samples("coh1", self.y_test1.index)
        age2 = self._get_age_series_for_samples("coh2", self.y_test2.index)

        y1 = (self.y_test1 == "Malignant").values.astype(int)
        y2 = (self.y_test2 == "Malignant").values.astype(int)

        for grp_name, mask1, mask2 in [
            ("Age_lt40", (age1 <= 40), (age2 <= 40)),
            ("Age_ge55", (age1 >= 55), (age2 >= 55)),
        ]:
            idx1 = np.where(mask1.values)[0]
            if len(idx1) > 0:
                m1 = self._calc_binary_metrics(y1[idx1], best_proba1[idx1], threshold=0.5)
                age_rows.append({
                     "Model": "SEM",          # ✅ 不用 BEST_SEM / name
                    "Subset": grp_name,
                    "Cohort": self.cohort_display["coh1"],
                    "n": int(len(idx1)),
                    "n_pos": int(y1[idx1].sum()),
                    **m1
                })

            idx2 = np.where(mask2.values)[0]
            if len(idx2) > 0:
                m2 = self._calc_binary_metrics(y2[idx2], best_proba2[idx2], threshold=0.5)
                age_rows.append({
                     "Model": "SEM",          # ✅ 不用 BEST_SEM / name
                    "Subset": grp_name,
                    "Cohort": self.cohort_display["coh2"],
                    "n": int(len(idx2)),
                    "n_pos": int(y2[idx2].sum()),
                    **m2
                })

        # 保存年龄子组结果
        if len(age_rows) > 0:
            age_df = pd.DataFrame(age_rows)
            age_path = os.path.join(self.output_dir, f"AgeSubgroup_metrics_{self.feature_mode}.csv")
            age_df.to_csv(age_path, index=False)
            print(f"✅ 年龄子组统计已保存: {age_path}")
        else:
            print("⚠️ 年龄子组统计：没有满足 Age<40 或 Age>=55 的样本，未生成文件。")

        print(f"\n✅ 多模型对比结果已保存到: {cmp_path}")

        # =========================
        # BMI 子组统计（SEM-only）
        # BMI <25 正常；BMI >=28 肥胖
        # =========================
        bmi_rows = []

        bmi1 = self._get_bmi_series_for_samples("coh1", self.y_test1.index)
        bmi2 = self._get_bmi_series_for_samples("coh2", self.y_test2.index)

        y1 = (self.y_test1 == "Malignant").values.astype(int)
        y2 = (self.y_test2 == "Malignant").values.astype(int)

        for grp_name, mask1, mask2 in [
            ("BMI_lt25_Normal", (bmi1 < 25), (bmi2 < 25)),
            ("BMI_ge28_Obese", (bmi1 >= 28), (bmi2 >= 28)),
        ]:
            # Cohort1(Test1)
            idx1 = np.where(mask1.values)[0]
            if len(idx1) > 0:
                m1 = self._calc_binary_metrics(y1[idx1], best_proba1[idx1], threshold=0.5)
                bmi_rows.append({
                    "Model": "SEM",
                    "Subset": grp_name,
                    "Cohort": self.cohort_display["coh1"],
                    "n": int(len(idx1)),
                    "n_pos": int(y1[idx1].sum()),
                    **m1
                })

            # Cohort2(External)
            idx2 = np.where(mask2.values)[0]
            if len(idx2) > 0:
                m2 = self._calc_binary_metrics(y2[idx2], best_proba2[idx2], threshold=0.5)
                bmi_rows.append({
                    "Model": "SEM",
                    "Subset": grp_name,
                    "Cohort": self.cohort_display["coh2"],
                    "n": int(len(idx2)),
                    "n_pos": int(y2[idx2].sum()),
                    **m2
                })

        # 保存 BMI 子组结果
        if len(bmi_rows) > 0:
            bmi_df = pd.DataFrame(bmi_rows)
            bmi_path = os.path.join(self.output_dir, f"BMISubgroup_metrics_{self.feature_mode}.csv")
            bmi_df.to_csv(bmi_path, index=False)
            print(f"✅ BMI 子组统计已保存: {bmi_path}")
        else:
            print("⚠️ BMI 子组统计：没有满足 BMI<25 或 BMI>=28 的样本，未生成文件。")

        return results_df

    def export_sem_final_metrics(self, threshold=0.5, n_bootstraps=1000, seed=40):
        """
        导出 SEM 最终模型在 Cohort1 (internal test) 和 Cohort2 (external) 上的
        完整分类指标：AUC / Accuracy / Precision / Recall / Specificity / F1 / NPV / MCC,
        每个指标含 95% bootstrap CI。

        前置条件：evaluate_models() 已经跑过 (self.additional_models["SEM"] 存在)。

        输出:
          {output_dir}/SEM_final_metrics_{feature_mode}.csv       (宽表, 一行一 cohort)
          {output_dir}/SEM_final_metrics_{feature_mode}_long.csv  (长表, 一行一 metric, 方便 ggplot)
        """
        from sklearn.metrics import (
            accuracy_score, precision_score, recall_score,
            f1_score, matthews_corrcoef, confusion_matrix
        )

        if "SEM" not in self.additional_models:
            print("⚠️ SEM 模型尚未训练，请先跑 evaluate_models()。")
            return None

        sem = self.additional_models["SEM"]

        cohorts = [
            ("Cohort1_Test1",    self.X_test1_scaled, self.y_test1),
            ("Cohort2_External", self.X_test2_scaled, self.y_test2),
        ]

        def _calc_one(yt, yp_proba, yp_pred):
            d = {}
            try:
                d["AUC"] = float(roc_auc_score(yt, yp_proba))
            except Exception:
                d["AUC"] = float("nan")
            d["Accuracy"]  = float(accuracy_score(yt, yp_pred))
            d["Precision"] = float(precision_score(yt, yp_pred, zero_division=0))
            d["Recall"]    = float(recall_score(yt, yp_pred, zero_division=0))
            d["F1"]        = float(f1_score(yt, yp_pred, zero_division=0))
            d["MCC"]       = float(matthews_corrcoef(yt, yp_pred))
            cm = confusion_matrix(yt, yp_pred, labels=[0, 1])
            if cm.shape == (2, 2):
                tn, fp, fn, tp = cm.ravel()
                d["Specificity"] = float(tn / (tn + fp)) if (tn + fp) > 0 else 0.0
                d["NPV"]         = float(tn / (tn + fn)) if (tn + fn) > 0 else 0.0
            else:
                d["Specificity"] = float("nan")
                d["NPV"]         = float("nan")
            return d

        metric_order = ["AUC", "Accuracy", "Precision", "Recall",
                        "Specificity", "F1", "NPV", "MCC"]

        long_rows = []
        wide_rows = []

        for cohort_name, X_scaled, y_series in cohorts:
            proba  = sem.predict_proba(X_scaled)[:, 1]
            y_bin  = (y_series == "Malignant").values.astype(int)
            y_pred = (proba >= threshold).astype(int)

            # 点估计
            point = _calc_one(y_bin, proba, y_pred)

            # Bootstrap CI
            rng = np.random.RandomState(seed)
            n = len(y_bin)
            boot = {k: [] for k in point.keys()}
            for _ in range(n_bootstraps):
                idx = rng.randint(0, n, n)
                if len(np.unique(y_bin[idx])) < 2:
                    continue
                try:
                    m = _calc_one(y_bin[idx], proba[idx], y_pred[idx])
                    for k, v in m.items():
                        if not np.isnan(v):
                            boot[k].append(v)
                except Exception:
                    continue

            wide_row = {
                "Feature_Mode": self.feature_mode,
                "Cohort": cohort_name,
                "N_total":    int(n),
                "N_positive": int(y_bin.sum()),
                "N_negative": int((y_bin == 0).sum()),
                "Threshold":  threshold,
                "n_bootstraps": int(n_bootstraps),
            }
            for k in metric_order:
                v  = point[k]
                lo = float(np.percentile(boot[k], 2.5))  if len(boot[k]) > 0 else float("nan")
                hi = float(np.percentile(boot[k], 97.5)) if len(boot[k]) > 0 else float("nan")
                wide_row[k]               = v
                wide_row[f"{k}_CI_low"]   = lo
                wide_row[f"{k}_CI_high"]  = hi
                long_rows.append({
                    "Feature_Mode": self.feature_mode,
                    "Cohort": cohort_name,
                    "Metric": k,
                    "Value":  v,
                    "CI_low": lo,
                    "CI_high": hi,
                })
            wide_rows.append(wide_row)

        wide_df = pd.DataFrame(wide_rows)
        long_df = pd.DataFrame(long_rows)

        wide_path = os.path.join(self.output_dir,
                                  f"SEM_final_metrics_{self.feature_mode}.csv")
        long_path = os.path.join(self.output_dir,
                                  f"SEM_final_metrics_{self.feature_mode}_long.csv")
        wide_df.to_csv(wide_path, index=False)
        long_df.to_csv(long_path, index=False)

        # 打印总结
        print(f"\n{'=' * 72}")
        print(f"SEM final metrics — feature_mode = {self.feature_mode}")
        print('=' * 72)
        for _, r in wide_df.iterrows():
            print(f"\n[{r['Cohort']}]  N = {int(r['N_total'])} "
                  f"({int(r['N_positive'])} Mal / {int(r['N_negative'])} Ben), "
                  f"threshold = {threshold}")
            for m in metric_order:
                print(f"  {m:12s} = {r[m]:.3f}  "
                      f"[{r[m+'_CI_low']:.3f}, {r[m+'_CI_high']:.3f}]")
        print(f"\n✅ Wide format: {wide_path}")
        print(f"✅ Long format: {long_path}")

        self.sem_final_metrics_wide = wide_df
        self.sem_final_metrics_long = long_df
        return wide_df

    def _export_predictions(self, model_name, cohort_tag, sample_index, y_series, proba, threshold=0.5):
        """
        导出逐样本预测：
        - sample_index: 样本名（Index）
        - y_series: 原始标签（Benign/Malignant）
        - proba: 预测概率（P(Malignant)）
        """
        y_true = (y_series == "Malignant").values.astype(int)
        y_pred = (proba >= threshold).astype(int)

        df = pd.DataFrame({
            "Sample": list(sample_index),
            "y_true": y_true,
            "y_label": list(y_series.values),
            "y_prob": proba,
            "y_pred": y_pred
        })

        out = os.path.join(self.output_dir, f"Pred_{model_name}_{cohort_tag}.csv")
        df.to_csv(out, index=False)
        print(f"✅ 已保存逐样本预测: {out}")
        return out

    def _get_age_series_for_samples(self, cohort_name, sample_index):
        """
        cohort_name: "coh1" / "coh2" / "coh3"
        返回与 sample_index 对齐的年龄 Series（float）
        """
        if cohort_name == "coh1":
            df = self.coh1_data
        elif cohort_name == "coh2":
            df = self.coh2_data
        elif cohort_name == "coh3":
            df = self.coh3_data
        else:
            raise ValueError("cohort_name must be one of ['coh1','coh2','coh3'].")

        age = df.loc[sample_index, "Age"].copy()
        age = pd.to_numeric(age, errors="coerce")
        return age

    def _get_bmi_series_for_samples(self, cohort_name, sample_index):
        """
        cohort_name: "coh1" / "coh2" / "coh3"
        返回与 sample_index 对齐的 BMI Series（float）
        """
        if cohort_name == "coh1":
            df = self.coh1_data
        elif cohort_name == "coh2":
            df = self.coh2_data
        elif cohort_name == "coh3":
            df = self.coh3_data
        else:
            raise ValueError("cohort_name must be one of ['coh1','coh2','coh3'].")

        bmi = df.loc[sample_index, "BMI"].copy()
        bmi = pd.to_numeric(bmi, errors="coerce")
        return bmi

    def _calc_binary_metrics(self, y_true, y_prob, threshold=0.5):
        """
        y_true: 0/1 numpy array
        y_prob: probability numpy array
        """
        y_pred = (y_prob >= threshold).astype(int)

        # 基础指标
        acc = accuracy_score(y_true, y_pred)
        prec = precision_score(y_true, y_pred, zero_division=0)
        rec = recall_score(y_true, y_pred, zero_division=0)
        f1v = f1_score(y_true, y_pred, zero_division=0)

        # confusion matrix / specificity / npv
        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan

        # AUC（需要两类）
        auc_val, ci_low, ci_high = np.nan, np.nan, np.nan
        if len(np.unique(y_true)) >= 2:
            auc_val = roc_auc_score(y_true, y_prob)
            ci_low, ci_high = bootstrap_auc(y_true, y_prob)

        return {
            "AUC": auc_val, "AUC_CI_low": ci_low, "AUC_CI_high": ci_high,
            "Accuracy": acc, "Precision": prec, "Recall": rec, "F1": f1v,
            "TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp),
            "Specificity": specificity, "NPV": npv
        }


    # ==================== SEM-only 评估（只验证 Cohort2/Test1） ====================
    def evaluate_sem_only_on_cohort2(self, save_tag=None, threshold=0.5):
        """
        仅训练 SEM，并在 Cohort2(=coh1 split出来的 Test1) 上评估。
        额外输出：Acc/Precision/Recall/F1/Specificity/NPV/Confusion matrix + AUC(含bootstrap CI)
        并保存性能文件到 output_dir。
        """
        # 确保已经做过 RFE + 标准化（feature_selection 会生成 X_train_scaled 等）
        if not hasattr(self, "X_train_scaled"):
            raise RuntimeError("请先运行 feature_selection(n_features=...) 以得到 X_train_scaled / X_test1_scaled。")

        # ---- 1) 定义 SEM（保持你原来的 SEM 配置不变） ----
        glm = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=3000, random_state=40,
        )
        lr = LogisticRegression(
            penalty="l2", solver="liblinear",
            max_iter=3000, random_state=40,
        )
        svm = SVC(kernel="linear", probability=True, random_state=40)
        gbm = GradientBoostingClassifier(
            learning_rate=0.01, n_estimators=200,
            max_depth=1, subsample=0.8, random_state=40
        )
        mlp = MLPClassifier(
            learning_rate_init=0.001, hidden_layer_sizes=(64, 32),
            activation="relu", alpha=50,
            solver="adam", max_iter=2000, random_state=40
        )

        base_estimators = [
            ("glm", glm),
            ("lr", lr),
            ("svm", svm),
            ("mlp", mlp),
            ("gbm", gbm),
        ]
        sem = StackingClassifier(
            estimators=base_estimators,
            final_estimator=LogisticRegression(
                penalty="l2", solver="liblinear",
                max_iter=3000, random_state=40,
            ),
            passthrough=False,
            n_jobs=1,
        )

        # ---- 2) 训练 SEM（训练集=coh1的train）----
        Xtr, Xva = self.X_train_scaled, self.X_test1_scaled
        ytr, yva = self.y_train, self.y_test1

        sem.fit(Xtr, ytr)
        self.additional_models["SEM"] = sem  # 仍然存着，方便你后续画图

        # ---- 3) 预测概率 & 二分类标签 ----
        proba = sem.predict_proba(Xva)[:, 1]
        y_true = (yva == "Malignant").values.astype(int)
        y_pred = (proba >= threshold).astype(int)

        # ---- 4) AUC + CI ----
        auc_val = roc_auc_score(y_true, proba)
        ci_low, ci_high = bootstrap_auc(y_true, proba)

        # ---- 5) 常用分类指标 ----
        acc = accuracy_score(y_true, y_pred)
        prec = precision_score(y_true, y_pred, zero_division=0)
        rec = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)

        tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
        specificity = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan  # negative predictive value

        # ---- 6) 汇总成一行性能表 ----
        row = {
            "Mode": self.feature_mode,
            "n_features": len(getattr(self, "selected_features", [])),
            "Cohort": self.cohort_display["coh1"],  # Cohort 2
            "Threshold": threshold,

            "AUC": auc_val,
            "AUC_CI_low": ci_low,
            "AUC_CI_high": ci_high,

            "Accuracy": acc,
            "Precision": prec,
            "Recall": rec,
            "F1": f1,

            "TN": int(tn), "FP": int(fp), "FN": int(fn), "TP": int(tp),
            "Specificity": specificity,
            "NPV": npv,
        }
        perf_df = pd.DataFrame([row])

        # ---- 7) 保存（保留性能文件）----
        tag = "" if save_tag is None else f"_{save_tag}"
        out_csv = os.path.join(self.output_dir, f"SEM_performance_cohort2{tag}.csv")
        perf_df.to_csv(out_csv, index=False)

        # 也把预测分数保存，方便你之后画别的曲线或做 DeLong 等
        pred_df = pd.DataFrame({
            "Sample": yva.index,
            "y_true": y_true,
            "y_prob": proba,
            "y_pred": y_pred,
        })
        out_pred = os.path.join(self.output_dir, f"SEM_predictions_cohort2{tag}.csv")
        pred_df.to_csv(out_pred, index=False)

        print(f"✅ SEM-only Cohort2 性能已保存: {out_csv}")
        print(f"✅ SEM-only Cohort2 预测已保存: {out_pred}")

        return perf_df



    # ==================== 4b. 不同特征数性能扫描 + 作图 ====================
    def sweep_feature_numbers_sem_only(
            self,
            feature_counts=(5, 10, 15, 20, 30, 40, 50),
            threshold=0.5
    ):
        """
        循环不同 RFE 特征数：
          - 每次：RFE选k个特征 -> 训练SEM -> 在 Cohort2(Test1) 验证
          - 每个k都单独保存性能文件
          - 最后汇总保存一份总表 + 画图用
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

            # 2) SEM-only 验证（Cohort2/Test1）
            perf_k = self.evaluate_sem_only_on_cohort2(save_tag=f"n{k}", threshold=threshold)
            all_rows.append(perf_k)

        sweep_df = pd.concat(all_rows, ignore_index=True)

        out_csv = os.path.join(self.output_dir, f"SEM_feature_number_sweep_{self.feature_mode}.csv")
        sweep_df.to_csv(out_csv, index=False)
        print(f"\n✅ SEM-only 不同特征数汇总结果已保存: {out_csv}")

        self.sem_feature_sweep_df = sweep_df
        return sweep_df

    def plot_sem_performance_vs_feature_number(self, sweep_df=None, ylim_auc=(0.75, 1.02)):
        """
        画 SEM 的性能随特征数变化：
          - 图1：AUC(含CI) vs 特征数
          - 图2：Acc/F1/Recall/Precision vs 特征数
        """
        if sweep_df is None:
            sweep_df = getattr(self, "sem_feature_sweep_df", None)
        if sweep_df is None:
            raise RuntimeError("请先运行 sweep_feature_numbers_sem_only() 得到 sweep_df。")

        sweep_df = sweep_df.sort_values("n_features")

        # ---- 图1：AUC vs 特征数（带CI）----
        fig1, ax1 = plt.subplots(1, 1, figsize=(7, 4))
        ax1.plot(sweep_df["n_features"], sweep_df["AUC"], marker="o", linewidth=2.6)
        ax1.fill_between(
            sweep_df["n_features"].values,
            sweep_df["AUC_CI_low"].values,
            sweep_df["AUC_CI_high"].values,
            alpha=0.15
        )
        ax1.set_xlabel("# Selected features (SVM-RFE)", fontweight="bold", fontsize=18)
        ax1.set_ylabel("AUC", fontweight="bold", fontsize=18)
        ax1.set_title(f"SEM performance vs #Features – {self.cohort_display['coh1']} ({self.feature_mode})",
                      fontweight="bold", fontsize=18)
        ax1.set_ylim(*ylim_auc)
        ax1.grid(True, alpha=0.3)
        plt.tight_layout()

        out1 = os.path.join(self.output_dir, f"SEM_AUC_vs_feature_number_{self.feature_mode}.png")
        fig1.savefig(out1, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()
        print(f"✅ 已保存: {out1}")

        # ---- 图2：Acc/F1/Recall/Precision vs 特征数 ----
        fig2, ax2 = plt.subplots(1, 1, figsize=(7, 5))
        for metric in ["Accuracy", "F1", "Recall", "Precision"]:
            if metric in sweep_df.columns:
                ax2.plot(sweep_df["n_features"], sweep_df[metric], marker="o", linewidth=2.3, label=metric)

        ax2.set_xlabel("# Selected features (SVM-RFE)", fontweight="bold", fontsize=14)
        ax2.set_ylabel("Score", fontweight="bold", fontsize=14)
        ax2.set_title(f"SEM classification metrics vs #Features – {self.cohort_display['coh1']} ({self.feature_mode})",
                      fontweight="bold", fontsize=14)
        ax2.set_ylim(0.7, 1.02)
        ax2.legend(loc="lower right")
        ax2.grid(True, alpha=0.3)
        plt.tight_layout()

        out2 = os.path.join(self.output_dir, f"SEM_metrics_vs_feature_number_{self.feature_mode}.png")
        fig2.savefig(out2, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()
        print(f"✅ 已保存: {out2}")

        return [fig1, fig2]



    def export_sem_roc_data(self):
        """
        导出 SEM 模型在 Cohort1 / Cohort2 / Early 子集上的 ROC 所需数据
        返回 dict，供跨模式（protein/metabolite/combined）汇总绘图
        """
        if "SEM" not in self.additional_models:
            # 确保已训练多模型（包含 SEM）
            self.evaluate_models()

        sem = self.additional_models["SEM"]

        Xtr, Xte1, Xte2 = self.X_train_scaled, self.X_test1_scaled, self.X_test2_scaled
        ytr, yte1, yte2 = self.y_train, self.y_test1, self.y_test2

        # Cohort1 / Cohort2 二分类标签
        y1_bin = (yte1 == "Malignant").values.astype(int)
        y2_bin = (yte2 == "Malignant").values.astype(int)

        # SEM 概率
        proba1 = sem.predict_proba(Xte1)[:, 1]
        proba2 = sem.predict_proba(Xte2)[:, 1]

        # ===== Early 子集（沿用你 plot_roc_curves 里的逻辑）=====
        idx_all_c1 = list(ytr.index) + list(yte1.index)
        y_all_c1 = pd.concat([ytr, yte1])
        early_mask_c1 = []
        early_labels_c1 = []

        for sample, label in zip(idx_all_c1, y_all_c1):
            stage = self.coh1_data.loc[sample, "Clinical stage"] if sample in self.coh1_data.index else np.nan
            if label == "Benign":
                early_mask_c1.append(True);
                early_labels_c1.append(0)
            elif (label == "Malignant" and pd.notna(stage) and stage in ["IA", "IB", "0", 0]):
                early_mask_c1.append(True);
                early_labels_c1.append(1)
            else:
                early_mask_c1.append(False);
                early_labels_c1.append(0)

        early_mask_c1 = np.array(early_mask_c1, dtype=bool)
        early_labels_c1 = np.array(early_labels_c1, dtype=int)

        proba_all_c1 = np.concatenate([
            sem.predict_proba(Xtr)[:, 1],
            sem.predict_proba(Xte1)[:, 1],
        ])
        y_early_c1 = early_labels_c1[early_mask_c1]
        proba_early_c1 = proba_all_c1[early_mask_c1]

        # Cohort2 early
        idx_all_c2 = list(yte2.index)
        y_all_c2 = yte2.copy()
        early_mask_c2 = []
        early_labels_c2 = []

        for sample, label in zip(idx_all_c2, y_all_c2):
            stage = self.coh2_data.loc[sample, "Clinical stage"] if sample in self.coh2_data.index else np.nan
            if label == "Benign":
                early_mask_c2.append(True);
                early_labels_c2.append(0)
            elif (label == "Malignant" and pd.notna(stage) and stage in ["IA", "IB", "0", 0]):
                early_mask_c2.append(True);
                early_labels_c2.append(1)
            else:
                early_mask_c2.append(False);
                early_labels_c2.append(0)

        early_mask_c2 = np.array(early_mask_c2, dtype=bool)
        early_labels_c2 = np.array(early_labels_c2, dtype=int)

        y_early_c2 = early_labels_c2[early_mask_c2]
        proba_early_c2 = proba2[early_mask_c2]

        return {
            "mode": self.feature_mode,
            "cohort1": {"y": y1_bin, "p": proba1},
            "cohort2": {"y": y2_bin, "p": proba2},
            "early_cohort1": {"y": y_early_c1, "p": proba_early_c1},
            "early_cohort2": {"y": y_early_c2, "p": proba_early_c2},
        }

    # ==================== 5. 多模型 ROC（逻辑基本沿用原来的） ====================
    def plot_roc_curves(self):
        """
        绘制ROC曲线：Cohort1 / Cohort2 多模型 + Early BC vs NBC
        """
        print("\n绘制ROC曲线...")
        display_mode_map = {
            "combined": "ProMeta-BC",
            "protein": "ProMeta-BC-Proteome",
            "metabolite": "ProMeta-BC-Metabolome"
        }

        if not getattr(self, "additional_models", None) or len(self.additional_models) == 0:
            print("尚未运行 evaluate_models()，自动调用该函数以获得多模型...")
            self.evaluate_models()

        preferred_order = ["GLM", "LR", "SVM", "GBM", "XGB", "MLP", "SEM"]
        model_names = [m for m in preferred_order if m in self.additional_models]

        color_list = [
            "#8B4CB8",
            "#2E8B57",
            "#FF6B35",
            "#2b93db",
            "#FFD700",
            "#00B0F0",
            "#7B68EE",
            "#FF1493",
        ]
        color_map = {name: color_list[i % len(color_list)] for i, name in enumerate(model_names)}

        Xtr, Xte1, Xte2 = self.X_train_scaled, self.X_test1_scaled, self.X_test2_scaled
        ytr, yte1, yte2 = self.y_train, self.y_test1, self.y_test2
        y1_bin = (yte1 == "Malignant").values
        y2_bin = (yte2 == "Malignant").values

        fig, axes = plt.subplots(1, 4, figsize=(24, 6))

        # 子图1：Cohort 1
        ax1 = axes[0]
        for name in model_names:
            clf = self.additional_models[name]
            proba1 = clf.predict_proba(Xte1)[:, 1]
            fpr1, tpr1, _ = roc_curve(y1_bin, proba1)
            auc1 = auc(fpr1, tpr1)
            ci_low1, ci_high1 = bootstrap_auc(y1_bin, proba1)

            ax1.plot(
                fpr1,
                tpr1,
                linewidth=2.5,
                color=color_map[name],
                label=f"{display_mode_map.get(name, name)}: AUC={auc1:.3f} (95% CI {ci_low1:.3f}-{ci_high1:.3f})",
            )

        ax1.plot([0, 1], [0, 1], "gray", linestyle="--", linewidth=1.5, alpha=0.8)
        ax1.set_xlim([0.0, 1.0])
        ax1.set_ylim([0.0, 1.05])
        ax1.set_xlabel("False positive rate", fontweight="bold")
        ax1.set_ylabel("True positive rate", fontweight="bold")
        ax1.set_title(f"Total malignant diagnosis - {self.cohort_display['coh1']} ({self.feature_mode})",
                      fontweight="bold", fontsize=14)
        ax1.legend(loc="lower right", fontsize=11)
        ax1.grid(True, alpha=0.3)

        # 子图2：Cohort 2
        ax2 = axes[1]
        for name in model_names:
            clf = self.additional_models[name]
            proba2 = clf.predict_proba(Xte2)[:, 1]
            fpr2, tpr2, _ = roc_curve(y2_bin, proba2)
            auc2 = auc(fpr2, tpr2)
            ci_low2, ci_high2 = bootstrap_auc(y2_bin, proba2)

            ax2.plot(
                fpr2,
                tpr2,
                linewidth=2.5,
                color=color_map[name],
                label=f"{display_mode_map.get(name, name)}: AUC={auc2:.3f} (95% CI {ci_low2:.3f}-{ci_high2:.3f})",
            )

        ax2.plot([0, 1], [0, 1], "gray", linestyle="--", linewidth=1.5, alpha=0.8)
        ax2.set_xlim([0.0, 1.0])
        ax2.set_ylim([0.0, 1.05])
        ax2.set_xlabel("False positive rate", fontweight="bold")
        ax2.set_ylabel("True positive rate", fontweight="bold")
        ax2.set_title(f"Total malignant diagnosis - {self.cohort_display['coh2']} ({self.feature_mode})",
                      fontweight="bold", fontsize=14)
        ax2.legend(loc="lower right", fontsize=11)
        ax2.grid(True, alpha=0.3)

        # ===== 早期 BC vs NBC 掩码构建，与原逻辑一致 =====
        idx_all_c1 = list(ytr.index) + list(yte1.index)
        y_all_c1 = pd.concat([ytr, yte1])
        early_mask_c1 = []
        early_labels_c1 = []

        for sample, label in zip(idx_all_c1, y_all_c1):
            if sample in self.coh1_data.index:
                stage = self.coh1_data.loc[sample, "Clinical stage"]
            else:
                stage = np.nan

            if label == "Benign":
                early_mask_c1.append(True)
                early_labels_c1.append(0)
            elif (
                    label == "Malignant"
                    and pd.notna(stage)
                    and stage in ["IA", "IB", "0", 0]
            ):
                early_mask_c1.append(True)
                early_labels_c1.append(1)
            else:
                early_mask_c1.append(False)
                early_labels_c1.append(0)

        early_mask_c1 = np.array(early_mask_c1, dtype=bool)
        early_labels_c1 = np.array(early_labels_c1, dtype=int)

        idx_all_c2 = list(yte2.index)
        y_all_c2 = yte2.copy()
        early_mask_c2 = []
        early_labels_c2 = []

        for sample, label in zip(idx_all_c2, y_all_c2):
            if sample in self.coh2_data.index:
                stage = self.coh2_data.loc[sample, "Clinical stage"]
            else:
                stage = np.nan

            if label == "Benign":
                early_mask_c2.append(True)
                early_labels_c2.append(0)
            elif (
                    label == "Malignant"
                    and pd.notna(stage)
                    and stage in ["IA", "IB", "0", 0]
            ):
                early_mask_c2.append(True)
                early_labels_c2.append(1)
            else:
                early_mask_c2.append(False)
                early_labels_c2.append(0)

        early_mask_c2 = np.array(early_mask_c2, dtype=bool)
        early_labels_c2 = np.array(early_labels_c2, dtype=int)

        X_all_c1 = np.vstack([Xtr, Xte1])

        # 子图3：Cohort 1 Early BC vs NBC
        ax3 = axes[2]
        if early_mask_c1.sum() > 0:
            print(
                f"Cohort 1 Early 分析：总有效样本 {early_mask_c1.sum()}，"
                f"其中 Early BC={early_labels_c1[early_mask_c1].sum()}，"
                f"NBC={early_mask_c1.sum() - early_labels_c1[early_mask_c1].sum()}"
            )

            for name in model_names:
                clf = self.additional_models[name]
                proba_all_c1 = np.concatenate([
                    clf.predict_proba(Xtr)[:, 1],
                    clf.predict_proba(Xte1)[:, 1],
                ])

                y_early_c1 = early_labels_c1[early_mask_c1]
                proba_early_c1 = proba_all_c1[early_mask_c1]

                if len(np.unique(y_early_c1)) < 2:
                    print(f"  模型 {name} 在 Cohort 1 Early 子集上只有单一类别，跳过 ROC。")
                    continue

                fpr3, tpr3, _ = roc_curve(y_early_c1, proba_early_c1)
                auc3 = auc(fpr3, tpr3)
                ci_low3, ci_high3 = bootstrap_auc(y_early_c1, proba_early_c1)

                ax3.plot(
                    fpr3,
                    tpr3,
                    linewidth=2.5,
                    color=color_map[name],
                    label=f"{display_mode_map.get(name, name)}: AUC={auc3:.3f} (95% CI {ci_low3:.3f}-{ci_high3:.3f})",
                )

            ax3.plot([0, 1], [0, 1], "gray", linestyle="--", linewidth=2, alpha=0.8)
            ax3.set_xlim([0.0, 1.0])
            ax3.set_ylim([0.0, 1.05])
            ax3.set_xlabel("False positive rate", fontweight="bold")
            ax3.set_ylabel("True positive rate", fontweight="bold")
            ax3.set_title(f"Early stage (0-I) diagnosis - {self.cohort_display['coh1']} ({self.feature_mode})",
                          fontweight="bold", fontsize=14)
            ax3.legend(loc="lower right", fontsize=12)
            ax3.grid(True, alpha=0.3)

        # 子图4：Cohort 2 Early BC vs NBC
        ax4 = axes[3]
        if early_mask_c2.sum() > 0:
            print(
                f"Cohort 2 Early 分析：总有效样本 {early_mask_c2.sum()}，"
                f"其中 Early BC={early_labels_c2[early_mask_c2].sum()}，"
                f"NBC={early_mask_c2.sum() - early_labels_c2[early_mask_c2].sum()}"
            )
            for name in model_names:
                clf = self.additional_models[name]
                proba_all_c2 = clf.predict_proba(Xte2)[:, 1]

                y_early_c2 = early_labels_c2[early_mask_c2]
                proba_early_c2 = proba_all_c2[early_mask_c2]

                if len(np.unique(y_early_c2)) < 2:
                    print(f"  模型 {name} 在 Cohort 2 Early 子集上只有单一类别，跳过 ROC。")
                    continue

                fpr4, tpr4, _ = roc_curve(y_early_c2, proba_early_c2)
                auc4 = auc(fpr4, tpr4)
                ci_low4, ci_high4 = bootstrap_auc(y_early_c2, proba_early_c2)

                ax4.plot(
                    fpr4,
                    tpr4,
                    linewidth=2.5,
                    color=color_map[name],
                    label=f"{display_mode_map.get(name, name)}: AUC={auc4:.3f} (95% CI {ci_low4:.3f}-{ci_high4:.3f})",
                )


            ax4.plot([0, 1], [0, 1], "gray", linestyle="--", linewidth=2, alpha=0.8)
            ax4.set_xlim([0.0, 1.0])
            ax4.set_ylim([0.0, 1.05])
            ax4.set_xlabel("False positive rate", fontweight="bold")
            ax4.set_ylabel("True positive rate", fontweight="bold")
            ax4.set_title(f"Early stage (0-I) diagnosis - {self.cohort_display['coh2']} ({self.feature_mode})",
                          fontweight="bold", fontsize=14)
            ax4.legend(loc="lower right", fontsize=11)
            ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        roc_path = f"{self.output_dir}/roc_curves_multimodel_with_early.png"
        fig.savefig(roc_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        print(f"多模型 ROC 曲线（含 Early BC）已保存到: {roc_path}")
        return fig

    # ==================== 6. 特征重要性 ====================
    # def plot_feature_importance(self):
    #     """绘制基于 SVM-RFE 的特征重要性（解释用）"""
    #     print("\n绘制特征重要性...")
    #
    #     fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    #
    #     top_features = self.feature_importance.head(10)
    #     colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(top_features)))
    #
    #     bars = ax.bar(
    #         range(len(top_features)),
    #         np.abs(top_features["coefficient"]),
    #         color=colors,
    #         alpha=0.8,
    #         edgecolor="black",
    #         linewidth=0.5,
    #     )
    #
    #     ax.set_xlabel("Features", fontsize=20, fontweight="bold")
    #     ax.set_ylabel("Contribution", fontsize=20, fontweight="bold")
    #     ax.set_title(f"Contribution for BC prediction ({self.feature_mode})",
    #                  fontsize=22, fontweight="bold", pad=20)
    #
    #     ax.set_xticks(range(len(top_features)))
    #
    #     def wrap_label(label, width=6):
    #         return textwrap.fill(label, width=width)
    #
    #     simplified_names = [
    #         wrap_label(
    #             f.split("_")[0] if "_" in f else f,
    #             width=10  # ✅ 每行最多 6 个字符，可改成 5 / 7
    #         )
    #         for f in top_features["feature"]
    #     ]
    #
    #     ax.tick_params(
    #         axis="both",
    #         which="major",
    #         labelsize=18,
    #         width=1.5,
    #         length=6
    #     )
    #     ax.set_xticklabels(simplified_names, rotation=45, ha="right", fontsize=16)
    #
    #     for bar, value in zip(bars, np.abs(top_features["coefficient"])):
    #         ax.text(
    #             bar.get_x() + bar.get_width() / 2.0,
    #             bar.get_height() + 0.001,
    #             f"{value:.3f}",
    #             ha="center",
    #             va="bottom",
    #             fontsize=18,
    #             fontweight="bold",
    #         )
    #     ax.grid(True, alpha=0.3, axis="y")
    #     ax.set_ylim(0, max(np.abs(top_features["coefficient"])) * 1.15)
    #
    #     plt.tight_layout()
    #
    #     importance_path = f"{self.output_dir}/feature_importance.png"
    #     fig.savefig(importance_path, dpi=300, bbox_inches="tight", facecolor="white")
    #     plt.show()
    #
    #     print(f"特征重要性图已保存到: {importance_path}")
    #     return fig

    def plot_feature_importance(self, top_n=None):
        """绘制基于 SVM-RFE 的特征重要性（默认展示全部；如需只看前N个可传 top_n=N）"""
        print("\n绘制特征重要性...")

        if not hasattr(self, "feature_importance") or self.feature_importance is None:
            raise RuntimeError("请先运行 feature_selection(n_features=...) 生成 feature_importance。")

        df = self.feature_importance.copy()

        # 默认展示全部（例如 n_features=20 就展示20个）
        if top_n is None:
            plot_df = df.copy()
        else:
            top_n = int(top_n)
            plot_df = df.head(top_n).copy()

        # 为了可视化：按绝对值从大到小
        plot_df["abs_coef"] = np.abs(plot_df["coefficient"])
        plot_df = plot_df.sort_values("abs_coef", ascending=True)  # barh 从下往上更美观

        # 自适应画布高度：每个特征0.55英寸，最低8英寸
        n = plot_df.shape[0]
        fig_h = max(8, 0.55 * n)
        fig, ax = plt.subplots(1, 1, figsize=(12, fig_h))

        # 颜色（viridis）
        colors = plt.cm.viridis(np.linspace(0.2, 0.85, n))

        # 标签 wrap
        def wrap_label(label, width=30):
            return textwrap.fill(str(label), width=width)

        y_labels = [wrap_label(f.split("_")[0] if "_" in f else f, width=30)
                    for f in plot_df["feature"].tolist()]

        y_pos = np.arange(n)

        bars = ax.barh(
            y_pos,
            plot_df["abs_coef"].values,
            color=colors,
            alpha=0.85,
            edgecolor="black",
            linewidth=0.6
        )

        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_labels, fontsize=28)
        for tick, feat in zip(ax.get_yticklabels(), plot_df["feature"].tolist()):
            if feat in self.keep_proteins:
                tick.set_color("#2b93db")
            elif feat in self.keep_metas:
                tick.set_color("#FF6B35")
            else:
                tick.set_color("black")
        ax.set_xlabel("Absolute coefficient (SVM-RFE)", fontsize=36, fontweight="bold")
        ax.set_title(f"Feature importance (n={n}) – {self.feature_mode}", fontsize=36, fontweight="bold", pad=12)

        ax.tick_params(axis="x", labelsize=28, width=1.5, length=6)
        ax.grid(True, alpha=0.25, axis="x")

        # 数值标注
        xmax = plot_df["abs_coef"].max() if n > 0 else 1.0
        for bar, val in zip(bars, plot_df["abs_coef"].values):
            ax.text(
                bar.get_width() + xmax * 0.01,
                bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}",
                va="center",
                ha="left",
                fontsize=28,
                fontweight="bold"
            )

        ax.set_xlim(0, xmax * 1.18)
        plt.tight_layout()

        # 文件名：如果 top_n=None 就表示全部（例如20）
        tag = f"all_{n}" if top_n is None else f"top_{top_n}"
        importance_path = os.path.join(self.output_dir, f"feature_importance_{tag}.png")
        fig.savefig(importance_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        print(f"特征重要性图已保存到: {importance_path}")
        return fig

    def plot_pca_analysis(self):
        """
        绘制 PCA + 预测值 分布图（按你原始逻辑）
        - 横轴：SEM 预测概率（Predicted value for BC）
        - 纵轴：PC1（PCA 从特征空间提取）
        PCA 仅在 Cohort2 的 Train+Test1 上拟合，再投影到 Cohort3（External）。
        """
        print("\n绘制PCA分析（x=SEM预测概率, y=PC1）...")

        # 确保已经训练多模型，并且包含 SEM
        if not getattr(self, "additional_models", None) or len(self.additional_models) == 0:
            print("尚未运行 evaluate_models()，自动调用该函数以获得多模型...")
            self.evaluate_models()

        if "SEM" not in self.additional_models:
            raise ValueError("当前 additional_models 中没有 'SEM' 模型，请检查 evaluate_models() 是否包含 SEM。")

        sem = self.additional_models["SEM"]

        # SEM 预测概率（横轴）
        proba_train_sem = sem.predict_proba(self.X_train_scaled)[:, 1]
        proba_test1_sem = sem.predict_proba(self.X_test1_scaled)[:, 1]
        proba_test2_sem = sem.predict_proba(self.X_test2_scaled)[:, 1]  # External / Cohort3

        # PCA：在 Train+Test1 拟合（避免使用 External 信息）
        X_all = np.vstack([self.X_train_scaled, self.X_test1_scaled])
        pca = PCA(n_components=2, random_state=40)
        X_pca = pca.fit_transform(X_all)

        X_train_pca = X_pca[:len(self.X_train_scaled)]
        X_test1_pca = X_pca[len(self.X_train_scaled):]
        X_test2_pca = pca.transform(self.X_test2_scaled)

        # 颜色
        benign_color = "#2E8B57"
        malignant_color = "#8B4CB8"

        cohort2_name = self.cohort_display["coh1"]  # Cohort 2
        cohort3_name = self.cohort_display["coh2"]  # Cohort 3

        legend_elements = [
            Patch(facecolor=benign_color, label="NBC"),
            Patch(facecolor=malignant_color, label="BC")
        ]

        # ----------------------------
        # d: Training set - by label
        # ----------------------------
        fig_d, ax_d = plt.subplots(1, 1, figsize=(8, 6))
        train_colors = [benign_color if lab == "Benign" else malignant_color for lab in self.y_train]

        ax_d.scatter(
            proba_train_sem, X_train_pca[:, 0],
            c=train_colors, alpha=0.75, s=60,
            edgecolors="black", linewidth=0.5
        )
        ax_d.axvline(x=0.5, color="gray", linestyle="--", alpha=0.7, linewidth=2)
        ax_d.set_xlabel("Predicted value for BC (SEM)", fontweight="bold", fontsize=20)
        ax_d.set_ylabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", fontweight="bold", fontsize=20)
        ax_d.set_title(f"Train+Validation set ({cohort2_name})", fontweight="bold", fontsize=22)
        ax_d.tick_params(
            axis="both",
            which="major",
            labelsize=18,
            width=1.5,
            length=6
        )
        ax_d.legend(handles=legend_elements, loc="upper right")
        ax_d.grid(True, alpha=0.3)

        plt.tight_layout()
        fig_d_path = f"{self.output_dir}/pca_train_labels.png"
        fig_d.savefig(fig_d_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        # ----------------------------
        # e: External (Cohort3) - by label
        # ----------------------------
        fig_e, ax_e = plt.subplots(1, 1, figsize=(8, 6))
        test2_colors = [benign_color if lab == "Benign" else malignant_color for lab in self.y_test2]

        ax_e.scatter(
            proba_test2_sem, X_test2_pca[:, 0],
            c=test2_colors, alpha=0.75, s=60,
            edgecolors="black", linewidth=0.5
        )
        ax_e.axvline(x=0.5, color="gray", linestyle="--", alpha=0.7, linewidth=2)
        ax_e.set_xlabel("Predicted value for BC (SEM)", fontweight="bold", fontsize=20)
        ax_e.set_ylabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", fontweight="bold", fontsize=20)
        ax_e.set_title(f"External set ({cohort3_name})", fontweight="bold", fontsize=24)

        ax_e.legend(handles=legend_elements, loc="upper right")
        ax_e.grid(True, alpha=0.3)
        ax_e.tick_params(
            axis="both",
            which="major",
            labelsize=18,
            width=1.5,
            length=6
        )
        plt.tight_layout()
        fig_e_path = f"{self.output_dir}/pca_external_labels.png"
        fig_e.savefig(fig_e_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        # ----------------------------
        # f1: Training set - by stage
        # ----------------------------
        fig_f1, ax_f1 = plt.subplots(1, 1, figsize=(8, 6))
        stage_colors_dict = {
            "NBC": benign_color,
            "Stage 0": "#FFB6C1",
            "Stage IA": "#FFA500",
            "Stage IB": "#FF8C00",
            "Stage II+": "#8B4513",
        }

        plotted_stages = set()
        for idx, sample in enumerate(self.y_train.index):
            if sample not in self.coh1_data.index:
                continue

            label = self.y_train.iloc[idx]
            stage_raw = self.coh1_data.loc[sample, "Clinical stage"]

            if label == "Benign":
                stage_label = "NBC"
            else:
                if pd.isna(stage_raw):
                    continue
                s = str(stage_raw).strip()
                if s in ["0"]:
                    stage_label = "Stage 0"
                elif s in ["IA"]:
                    stage_label = "Stage IA"
                elif s in ["IB"]:
                    stage_label = "Stage IB"
                else:
                    stage_label = "Stage II+"

            color = stage_colors_dict.get(stage_label, "#777777")
            show_label = stage_label if stage_label not in plotted_stages else ""
            plotted_stages.add(stage_label)

            ax_f1.scatter(
                proba_train_sem[idx], X_train_pca[idx, 0],
                c=color, alpha=0.75, s=60,
                label=show_label,
                edgecolors="black", linewidth=0.5
            )

        ax_f1.axvline(x=0.5, color="gray", linestyle="--", alpha=0.7, linewidth=2)
        ax_f1.set_xlabel("Predicted value for BC (SEM)", fontweight="bold", fontsize=20)
        ax_f1.set_ylabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", fontweight="bold", fontsize=20)
        ax_f1.set_title(f"Train+Validation set - by stage ({cohort2_name})", fontweight="bold", fontsize=22)
        ax_f1.tick_params(
            axis="both",
            which="major",
            labelsize=18,
            width=1.5,
            length=6
        )
        ax_f1.legend(loc="upper right")
        ax_f1.grid(True, alpha=0.3)

        plt.tight_layout()
        fig_f1_path = f"{self.output_dir}/pca_train_stages.png"
        fig_f1.savefig(fig_f1_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        # ----------------------------
        # f2: External (Cohort3) - by stage
        # ----------------------------
        fig_f2, ax_f2 = plt.subplots(1, 1, figsize=(8, 6))

        plotted_stages = set()
        for idx, sample in enumerate(self.y_test2.index):
            if sample not in self.coh2_data.index:
                continue

            label = self.y_test2.iloc[idx]
            stage_raw = self.coh2_data.loc[sample, "Clinical stage"]

            if label == "Benign":
                stage_label = "NBC"
            else:
                if pd.isna(stage_raw):
                    continue
                s = str(stage_raw).strip()
                if s in ["0"]:
                    stage_label = "Stage 0"
                elif s in ["IA"]:
                    stage_label = "Stage IA"
                elif s in ["IB"]:
                    stage_label = "Stage IB"
                else:
                    stage_label = "Stage II+"

            color = stage_colors_dict.get(stage_label, "#777777")
            show_label = stage_label if stage_label not in plotted_stages else ""
            plotted_stages.add(stage_label)

            ax_f2.scatter(
                proba_test2_sem[idx], X_test2_pca[idx, 0],
                c=color, alpha=0.75, s=60,
                label=show_label,
                edgecolors="black", linewidth=0.5
            )

        ax_f2.axvline(x=0.5, color="gray", linestyle="--", alpha=0.7, linewidth=2)
        ax_f2.set_xlabel("Predicted value for BC (SEM)", fontweight="bold", fontsize=20)
        ax_f2.set_ylabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", fontweight="bold", fontsize=20)
        ax_f2.set_title(f"External set - by stage ({cohort3_name})", fontweight="bold", fontsize=22)
        ax_f2.tick_params(
            axis="both",
            which="major",
            labelsize=18,
            width=1.5,
            length=6
        )
        ax_f2.legend(loc="upper right")
        ax_f2.grid(True, alpha=0.3)

        plt.tight_layout()
        fig_f2_path = f"{self.output_dir}/pca_external_stages.png"
        fig_f2.savefig(fig_f2_path, dpi=300, bbox_inches="tight", facecolor="white")
        plt.show()

        print("✅ PCA + Predicted value 图已保存：")
        print(f"- {fig_d_path}")
        print(f"- {fig_e_path}")
        print(f"- {fig_f1_path}")
        print(f"- {fig_f2_path}")

        return [fig_d, fig_e, fig_f1, fig_f2]

    # ==================== 8. 主流程 ====================
    def run_complete_analysis(self):
        """运行完整分析（本模式）"""
        print(f"\n开始完整的乳腺癌特征分析 - 模式: {self.feature_mode}")
        print("=" * 60)

        # 1. 数据加载
        self.load_data()

        # 2. 特征准备
        self.prepare_features()

        # 3. (新增) 训练集 5-fold CV 选方法（可选，但你现在想要的就是这个）
        self.select_model_by_cv_on_training(n_splits=5, n_features=20)

        # 4. 然后再用全 training set 做一次 RFE，固定面板
        self.feature_selection(n_features=20)

        # 5. 再 evaluate_models() 去 internal test + external 报告性能
        perf_df = self.evaluate_models()

        # 5.1 导出 SEM final model 的完整 metrics (AUC/F1/ACC/Precision/Recall/Spec/NPV/MCC + 95% CI)
        self.export_sem_final_metrics()

        # 5. 可视化
        self.plot_roc_curves()
        # self.plot_feature_importance()
        self.plot_feature_importance(top_n=None)  # None=全部(即20个)
        self.plot_pca_analysis()

        print("\n=== 各模型性能列表 ===")
        print(perf_df)

        print("\n=== 分析完成 ===")
        print(f"所有结果已保存到: {self.output_dir}/")

        return {
            'selected_features': self.selected_features,
            'feature_importance': self.feature_importance,
            'model_performance': perf_df,
            'sem_roc_data': self.export_sem_roc_data(),  # ✅新增
        }



def plot_sem_mode_roc_comparison(results_all, outdir="./breast_cancer_figures_SEM_compare"):
    os.makedirs(outdir, exist_ok=True)

    # 取三种模式的 SEM ROC 数据
    sem_data = {m: results_all[m]["sem_roc_data"] for m in ["combined", "protein", "metabolite"]}

    # 颜色
    colors = {
        "combined": "#8B4CB8",
        "protein": "#2b93db",
        "metabolite": "#FF6B35"
    }

    # (key, title) 按 2×2 顺序排
    panels = [
        ("cohort1", "Total malignant diagnosis – Cohort 2"),
        ("cohort2", "Total malignant diagnosis – Cohort 3"),
        ("early_cohort1", "Early stage (0-I) diagnosis – Cohort 2"),
        ("early_cohort2", "Early stage (0-I) diagnosis – Cohort 3"),
    ]
    mode_display = {
        "combined": "ProMeta-BC",
        "protein": "ProMeta-BC-Proteome",
        "metabolite": "ProMeta-BC-Metabolome"
    }

    # ✅ 2 × 2 布局
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()  # 变成一维，方便循环

    for ax, (key, title) in zip(axes, panels):
        for mode in ["combined", "protein", "metabolite"]:
            y = sem_data[mode][key]["y"]
            p = sem_data[mode][key]["p"]

            # Early 子集可能只有单一类别
            if len(np.unique(y)) < 2:
                continue

            fpr, tpr, _ = roc_curve(y, p)
            this_auc = auc(fpr, tpr)
            ci_low, ci_high = bootstrap_auc(y, p)

            ax.plot(
                fpr, tpr,
                linewidth=2.8,
                color=colors[mode],
                label=f"{mode_display[mode]}:\nAUC={this_auc:.3f} (95% CI {ci_low:.3f}-{ci_high:.3f})"
            )

        ax.plot([0, 1], [0, 1], "gray", linestyle="--", linewidth=1.6, alpha=0.8)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("False positive rate", fontweight="bold", fontsize=18)
        ax.set_ylabel("True positive rate", fontweight="bold", fontsize=18)
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=14,
            width=1.5,
            length=6
        )
        ax.set_title(title, fontweight="bold", fontsize=20)
        ax.legend(loc="lower right", fontsize=16)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    save_path = os.path.join(outdir, "SEM_only_ROC_compare_modes_2x2.png")
    fig.savefig(save_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()

    print(f"✅ SEM-only 三模式 ROC（2×2）对比图已保存: {save_path}")
    return fig

def export_all_models_auc_long_format(
    results_all,
    out_csv="./breast_cancer_figures_SEM_compare/all_models_auc_long_format.csv",
    cohort_label="Internal test (Cohort 2 split, Test1)"
):
    """
    仅汇总三种 feature_mode 下所有模型在 Cohort1(Test1, internal test) 的 AUC
    输出 long format:
      Model | AUC | Feature_Mode | Cohort
    """
    rows = []
    for mode in ["combined", "protein", "metabolite"]:
        perf_df = results_all[mode]["model_performance"].copy()

        need_cols = {"Model", "Cohort1_AUC"}
        if not need_cols.issubset(set(perf_df.columns)):
            raise RuntimeError(f"model_performance 缺少必要列: {need_cols - set(perf_df.columns)}")

        # 只保留主模型对比行：有 Cohort1_AUC 的那些
        perf_main = perf_df.loc[:, ["Model", "Cohort1_AUC"]].dropna(subset=["Cohort1_AUC"])
        perf_main = perf_main.drop_duplicates(subset=["Model"], keep="first")

        for _, r in perf_main.iterrows():
            rows.append({
                "Model": r["Model"],
                "AUC": float(r["Cohort1_AUC"]),
                "Feature_Mode": mode,
                "Cohort": cohort_label
            })

    out_df = pd.DataFrame(rows, columns=["Model", "AUC", "Feature_Mode", "Cohort"])
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    out_df.to_csv(out_csv, index=False)
    print(f"✅ 已保存（仅 internal test）all models AUC long format: {out_csv}")
    return out_df


def main():
    """
    一次性跑三套模型：
    - combined: 蛋白 + 代谢物
    - protein: 仅蛋白
    - metabolite: 仅代谢物
    每套结果单独输出到 ./breast_cancer_figures_<mode> 下
    """
    results_all = {}

    for mode in ["combined", "protein", "metabolite"]:
        outdir = f"./breast_cancer_figures_{mode}"
        analyzer = BreastCancerAnalysis(output_dir=outdir, feature_mode=mode)
        results_all[mode] = analyzer.run_complete_analysis()

        # ✅ 新增：仅 SEM，在 Cohort2(Test1) 上做不同特征数验证 + 保存性能文件 + 作图
        sweep_df = analyzer.sweep_feature_numbers_sem_only(
            feature_counts=range(5, 35),
            threshold=0.5
        )
        analyzer.plot_sem_performance_vs_feature_number(sweep_df=sweep_df)

    # ✅ 新增：只用 SEM，画三种模式（protein/metabolite/combined）的 ROC 对比
    plot_sem_mode_roc_comparison(results_all, outdir="./breast_cancer_figures_SEM_compare")

    # ✅ 新增：导出所有模型 AUC long format（与你给的csv格式一致）
    export_all_models_auc_long_format(
        results_all,
        out_csv="./breast_cancer_figures_SEM_compare/all_models_auc_long_format.csv",
        cohort_label = "Internal test (Cohort 2 split, Test1)"
    )

    return results_all


if __name__ == "__main__":
    results = main()
