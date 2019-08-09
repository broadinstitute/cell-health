import os
import sys
import pytest
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNet, SGDClassifier
from sklearn.pipeline import Pipeline

abs_path = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(abs_path, ".."))
from scripts.ml_utils import load_train_test, CellHealthPredict

np.random.seed(123)

# Define mock data
x_df = pd.DataFrame(
    [
        [x * 0.5 for x in [1, 5, 7, 2, 3, 9, 1]],
        [x * 0.1 for x in [1, 5, 7, 2, 3, 9, 1]],
        [x * 0.2 for x in [1, 5, 7, 2, 3, 9, 1]],
        [1, 5, 7, 2, 3, 9, 1],
        [x * 0.1 for x in [7, 0, 1, 3, 1, 0, 5]],
        [x * 0.2 for x in [7, 0, 1, 3, 1, 0, 5]],
        [x * 0.7 for x in [7, 0, 1, 3, 1, 0, 5]],
        [7, 0, 1, 3, 1, 0, 5],
    ],
    index=["a", "b", "c", "d", "e", "f", "g", "h"],
    columns=["x", "y", "z", "xx", "yy", "zz", "p"],
)

y_df = pd.DataFrame(
    [
        [x * 0.5 for x in [1, 1, 1, 0.01, 0.01, 0.01, 0.01]],
        [x * 1.1 for x in [1, 1, 1, 0.01, 0.01, 0.01, 0.01]],
        [x * 0.2 for x in [1, 1, 1, 0.01, 0.01, 0.01, 0.01]],
        [1, 1, 1, 0.01, 0.01, 0.01, 0.01],
        [x * 0.5 for x in [0.01, 0.01, 0.01, 1, 1, 1, 1]],
        [x * 1.2 for x in [0.01, 0.01, 0.01, 1, 1, 1, 1]],
        [x * 0.7 for x in [0.01, 0.01, 0.01, 1, 1, 1, 1]],
        [0.01, 0.01, 0.01, 1, 1, 1, 1],
    ],
    index=["a", "b", "c", "d", "e", "f", "g", "h"],
    columns=["target_{}".format(x) for x in ["x", "y", "z", "xx", "yy", "zz", "p"]],
)
y_df.loc["d", "target_p"] = np.nan

alphas = [0.01]
l1_ratios = [0.1]
n_folds = 2

regression_parameters = {"regress__alpha": alphas, "regress__l1_ratio": l1_ratios}

clf_parameters = {
    "classify__loss": ["log"],
    "classify__penalty": ["elasticnet"],
    "classify__alpha": alphas,
    "classify__l1_ratio": l1_ratios,
}

estimator_regressor = Pipeline(
    steps=[("regress", ElasticNet(random_state=42, max_iter=2000, tol=1e-3))]
)

estimator_classifier = Pipeline(
    steps=[
        (
            "classify",
            SGDClassifier(
                random_state=0,
                class_weight="balanced",
                max_iter=2000,
                shuffle=True,
                tol=1e-3,
            ),
        )
    ]
)

chp = CellHealthPredict(
    x_df=x_df,
    y_df=y_df,
    parameters=clf_parameters,
    estimator=estimator_classifier,
    n_folds=n_folds,
    cv_scoring="roc_auc",
)

chp_reg = CellHealthPredict(
    x_df=x_df,
    y_df=y_df,
    parameters=regression_parameters,
    estimator=estimator_regressor,
    n_folds=n_folds,
    cv_scoring="r2",
)


class TestCellHealthPredictClassify(object):
    def test_attr(self):
        assert hasattr(chp, "x_df")
        assert hasattr(chp, "y_df")
        assert hasattr(chp, "parameters")
        assert hasattr(chp, "estimator")
        assert hasattr(chp, "n_folds")
        assert hasattr(chp, "cv_scoring")
        assert hasattr(chp, "shuffle_key")
        assert hasattr(chp, "named_step")
        assert hasattr(chp, "return_train_score")
        assert hasattr(chp, "profile_ids")
        assert hasattr(chp, "is_fit")
        assert hasattr(chp, "kmeans_fit")
        assert hasattr(chp, "cv_pipeline")

    def test_fit_notarget(self):
        with pytest.raises(AssertionError) as e:
            chp.fit_cell_health_target("nothing")

        assert str(e.value) == "Feature must exist in y_df"

    def test_set_target(self):
        chp.set_target(target="target_x")
        assert hasattr(chp, "target")

    def test_set_y_transform(self):
        chp.set_y_transform(y_transform="binarize")
        assert hasattr(chp, "y_transform")

    def test_raw_recode(self):
        chp.set_target(target="target_zz")
        chp.set_y_transform(y_transform="raw")
        y = chp.recode_y(y=chp.y_df)

        y_expect = pd.Series(
            [0.005, 0.011, 0.002, 0.01, 0.5, 1.2, 0.7, 1], name="target_zz", dtype=float
        )
        y_expect.index = ["a", "b", "c", "d", "e", "f", "g", "h"]

        pd.testing.assert_series_equal(y, y_expect, check_dtype=False)

    def test_binarize_recode(self):
        y = chp.recode_y(y=chp.y_df, target="target_zz", y_transform="binarize")

        y_expect = pd.Series([0, 0, 0, 0, 1, 1, 1, 1], name="target_zz", dtype=np.int32)
        y_expect.index = ["a", "b", "c", "d", "e", "f", "g", "h"]

        pd.testing.assert_series_equal(y, y_expect)

    def test_zeroone_recode(self):
        y = chp.recode_y(y=chp.y_df, target="target_zz", y_transform="zero-one")

        y_expect = pd.Series(
            [0, 0.01, 0.0, 0.01, 0.42, 1, 0.58, 0.83], name="target_zz", dtype=float
        )
        y_expect.index = ["a", "b", "c", "d", "e", "f", "g", "h"]

        pd.testing.assert_series_equal(y.round(2), y_expect)

    def test_realign_missing_data(self):
        x_r, y_r, n = chp.realign_missing_data(chp.x_df, chp.y_df.loc[:, "target_p"])

        x_r_expect = pd.DataFrame(
            [
                [x * 0.5 for x in [1, 5, 7, 2, 3, 9, 1]],
                [x * 0.1 for x in [1, 5, 7, 2, 3, 9, 1]],
                [x * 0.2 for x in [1, 5, 7, 2, 3, 9, 1]],
                [x * 0.1 for x in [7, 0, 1, 3, 1, 0, 5]],
                [x * 0.2 for x in [7, 0, 1, 3, 1, 0, 5]],
                [x * 0.7 for x in [7, 0, 1, 3, 1, 0, 5]],
                [7, 0, 1, 3, 1, 0, 5],
            ],
            index=["a", "b", "c", "e", "f", "g", "h"],
            columns=["x", "y", "z", "xx", "yy", "zz", "p"],
        )
        pd.testing.assert_frame_equal(x_r, x_r_expect)

        y_r_expect = pd.Series(
            [0.005, 0.011, 0.002, 0.5, 1.2, 0.7, 1], name="target_p", dtype=float
        )
        y_r_expect.index = ["a", "b", "c", "e", "f", "g", "h"]
        pd.testing.assert_series_equal(y_r, y_r_expect)

        assert n == 1

    def test_fit_cell_health_target(self):
        output = chp.fit_cell_health_target(target="target_zz", y_transform="binarize")

        assert chp.is_fit
        assert output

    def test_get_cv_results(self):
        output = chp.fit_cell_health_target(target="target_zz", y_transform="binarize")

        cv = chp.get_cv_results()

        cv_expect = pd.DataFrame(
            [0.1, 0.01, 1.0, "target_zz", "shuffle_false", "binarize"]
        ).transpose()

        cv_expect.columns = [
            "l1_ratio",
            "alpha",
            "mean_test_score",
            "target",
            "shuffle",
            "y_transform",
        ]

        pd.testing.assert_frame_equal(cv, cv_expect, check_dtype=False)

    def test_binarize_predict(self):
        output = chp.fit_cell_health_target(target="target_zz", y_transform="binarize")

        pred = chp.predict(x_df.iloc[range(2, 5), :])
        pred_expect = pd.Series([0, 0, 1], name="target_zz", dtype=np.int32)

        pd.testing.assert_series_equal(pred, pred_expect)

    def test_binarize_get_coefficients(self):
        output = chp.fit_cell_health_target(target="target_zz", y_transform="binarize")

        coef = chp.get_coefficients().round(2)

        coef_exp = pd.DataFrame(
            [
                ["x", 7.75, 7.75, "target_zz", "binarize", "shuffle_false"],
                ["p", 5.43, 5.43, "target_zz", "binarize", "shuffle_false"],
                ["xx", 2.80, 2.80, "target_zz", "binarize", "shuffle_false"],
                ["zz", -2.70, 2.70, "target_zz", "binarize", "shuffle_false"],
                ["y", -1.47, 1.47, "target_zz", "binarize", "shuffle_false"],
                ["z", -0.92, 0.92, "target_zz", "binarize", "shuffle_false"],
                ["yy", 0.17, 0.17, "target_zz", "binarize", "shuffle_false"],
            ]
        )
        coef_exp.columns = [
            "feature",
            "weight",
            "abs_weight",
            "target",
            "y_transform",
            "shuffle",
        ]

        pd.testing.assert_frame_equal(coef, coef_exp)

    def test_binarize_get_performance(self):
        output = chp.fit_cell_health_target(target="target_zz", y_transform="binarize")

        roc_df, pr_df, y_true, y_pred = chp.get_performance(return_y=True)

        roc_df_expect = pd.DataFrame(
            [
                [
                    0.0,
                    0.0,
                    2,
                    "roc",
                    "target_zz",
                    1.0,
                    "train",
                    "shuffle_false",
                    "binarize",
                    4,
                ],
                [
                    0.0,
                    1.0,
                    1,
                    "roc",
                    "target_zz",
                    1.0,
                    "train",
                    "shuffle_false",
                    "binarize",
                    4,
                ],
                [
                    1.0,
                    1.0,
                    0,
                    "roc",
                    "target_zz",
                    1.0,
                    "train",
                    "shuffle_false",
                    "binarize",
                    4,
                ],
            ]
        )
        roc_df_expect.columns = [
            "fpr",
            "tpr",
            "threshold",
            "metric",
            "target",
            "auc",
            "data_fit",
            "shuffle",
            "y_transform",
            "min_class_count",
        ]

        pd.testing.assert_frame_equal(roc_df, roc_df_expect, check_dtype=False)

        pr_df_expect = pd.DataFrame(
            [
                [
                    1.0,
                    1.0,
                    1.0,
                    "aupr",
                    "target_zz",
                    1.0,
                    "train",
                    "shuffle_false",
                    "binarize",
                    4,
                ],
                [
                    1.0,
                    0.0,
                    np.nan,
                    "aupr",
                    "target_zz",
                    1.0,
                    "train",
                    "shuffle_false",
                    "binarize",
                    4,
                ],
            ]
        )
        pr_df_expect.columns = [
            "precision",
            "recall",
            "threshold",
            "metric",
            "target",
            "auc",
            "data_fit",
            "shuffle",
            "y_transform",
            "min_class_count",
        ]

        pd.testing.assert_frame_equal(pr_df, pr_df_expect)

        y_expect = pd.DataFrame(
            [["a", "b", "c", "d", "e", "f", "g", "h"], [0, 0, 0, 0, 1, 1, 1, 1]]
        ).transpose()
        y_expect = y_expect.assign(
            target="target_zz",
            data_type="train",
            shuffle="shuffle_false",
            y_transform="binarize",
        )
        y_expect.columns = [
            "Metadata_profile_id",
            "recode_target_value",
            "target",
            "data_type",
            "shuffle",
            "y_transform",
        ]
        pd.testing.assert_frame_equal(
            y_expect.assign(y_type="y_true"), y_true, check_dtype=False
        )
        pd.testing.assert_frame_equal(
            y_expect.assign(y_type="y_pred"), y_pred, check_dtype=False
        )


class TestCellHealthPredictRegression(object):
    def test_attr(self):
        assert hasattr(chp_reg, "x_df")
        assert hasattr(chp_reg, "y_df")
        assert hasattr(chp_reg, "parameters")
        assert hasattr(chp_reg, "estimator")
        assert hasattr(chp_reg, "n_folds")
        assert hasattr(chp_reg, "cv_scoring")
        assert hasattr(chp_reg, "shuffle_key")
        assert hasattr(chp_reg, "named_step")
        assert hasattr(chp_reg, "return_train_score")
        assert hasattr(chp_reg, "profile_ids")
        assert hasattr(chp_reg, "is_fit")
        assert hasattr(chp_reg, "kmeans_fit")
        assert hasattr(chp_reg, "cv_pipeline")

    def test_fit_cell_health_target(self):
        output = chp_reg.fit_cell_health_target(
            target="target_zz", y_transform="binarize"
        )

        assert chp_reg.is_fit
        assert output

    def test_binarize_get_coefficients(self):
        output = chp_reg.fit_cell_health_target(target="target_zz", y_transform="raw")

        coef = chp_reg.get_coefficients().round(2)

        coef_exp = pd.DataFrame(
            [
                ["x", 0.09, 0.09, "target_zz", "raw", "shuffle_false"],
                ["zz", -0.04, 0.04, "target_zz", "raw", "shuffle_false"],
                ["z", -0.03, 0.03, "target_zz", "raw", "shuffle_false"],
                ["p", 0.02, 0.02, "target_zz", "raw", "shuffle_false"],
                ["y", -0.00, 0.00, "target_zz", "raw", "shuffle_false"],
                ["xx", 0.00, 0.00, "target_zz", "raw", "shuffle_false"],
                ["yy", -0.00, 0.00, "target_zz", "raw", "shuffle_false"],
            ]
        )
        coef_exp.columns = [
            "feature",
            "weight",
            "abs_weight",
            "target",
            "y_transform",
            "shuffle",
        ]

        pd.testing.assert_frame_equal(coef, coef_exp)

    def test_raw_recode(self):
        y = chp_reg.recode_y(y=chp_reg.y_df, target="target_zz", y_transform="raw")

        pd.testing.assert_series_equal(y, chp_reg.y_df.target_zz)

    def test_reg_get_performance(self):
        output = chp_reg.fit_cell_health_target(target="target_p", y_transform="raw")

        mse_df, r2_df, y_true, y_pred = chp_reg.get_performance(return_y=True)

        mse_df_expect = pd.DataFrame(
            [0.08, "mse", "target_p", "train", "shuffle_false", "raw"]
        ).transpose()
        mse_df_expect.columns = [
            "mse",
            "metric",
            "target",
            "data_fit",
            "shuffle",
            "y_transform",
        ]

        pd.testing.assert_frame_equal(mse_df.round(2), mse_df_expect, check_dtype=False)

        r2_df_expect = pd.DataFrame(
            [0.61, "r_two", "target_p", "train", "shuffle_false", "raw"]
        ).transpose()
        r2_df_expect.columns = [
            "mse",
            "metric",
            "target",
            "data_fit",
            "shuffle",
            "y_transform",
        ]

        pd.testing.assert_frame_equal(r2_df.round(2), r2_df_expect, check_dtype=False)

        y_expect = pd.DataFrame(
            [
                ["a", "b", "c", "e", "f", "g", "h"],
                [0.005, 0.011, 0.002, 0.5, 1.2, 0.7, 1.0],
            ]
        ).transpose()
        y_expect = y_expect.assign(
            target="target_p",
            data_type="train",
            shuffle="shuffle_false",
            y_transform="raw",
        )
        y_expect.columns = [
            "Metadata_profile_id",
            "recode_target_value",
            "target",
            "data_type",
            "shuffle",
            "y_transform",
        ]
        pd.testing.assert_frame_equal(
            y_expect.assign(y_type="y_true"), y_true, check_dtype=False
        )

        y_expect.recode_target_value = [-0.15, 0.36, 0.23, 0.53, 0.59, 0.85, 1.01]
        pd.testing.assert_frame_equal(
            y_expect.assign(y_type="y_pred"), y_pred.round(2), check_dtype=False
        )

    def test_reg_get_performance_test_set(self):
        output = chp_reg.fit_cell_health_target(target="target_p", y_transform="raw")

        mse_df, r2_df, y_true, y_pred = chp_reg.get_performance(return_y=True)
        mse_test_df, r2_test_df, y_test_true, y_test_pred = chp_reg.get_performance(
            return_y=True, x_test=x_df, y_test=y_df
        )

        # Recode "train" to "test"
        mse_df.data_fit = "test"
        r_two_df.data_fit = "test"
        y_true.data_type = "test"
        y_pred.data_type = "test"

        pd.testing.assert_frame_equal(mse_df, mse_test_df)
        pd.testing.assert_frame_equal(r2_df, r2_test_df)
        pd.testing.assert_frame_equal(y_true, y_test_true)
        pd.testing.assert_frame_equal(y_pred, y_test_pred)
