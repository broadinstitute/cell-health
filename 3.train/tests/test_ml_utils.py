import os
import sys
import pytest
import numpy as np
import pandas as pd
from sklearn.linear_model import SGDRegressor, SGDClassifier
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

regression_parameters = {
    "regress__loss": ["squared_loss"],
    "regress__penalty": ["elasticnet"],
    "regress__alpha": alphas,
    "regress__l1_ratio": l1_ratios,
}

clf_parameters = {
    "classify__loss": ["log"],
    "classify__penalty": ["elasticnet"],
    "classify__alpha": alphas,
    "classify__l1_ratio": l1_ratios,
}
estimator_regressor = Pipeline(
    steps=[("regress", SGDRegressor(max_iter=2000, random_state=42, tol=1e-3))]
)

estimator_classifier = Pipeline(
    steps=[
        (
            "classify",
            SGDClassifier(
                random_state=0, class_weight="balanced", max_iter=2000, tol=1e-3
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


class TestCellHealthPredict(object):
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
        y = chp.recode_y(y=chp.y_df, target="target_zz", y_transform="raw")

        y_expect = pd.Series(
            [0.005, 0.011, 0.002, 0.01, 0.5, 1.2, 0.7, 1], name="target_zz", dtype=float
        )
        y_expect.index = ["a", "b", "c", "d", "e", "f", "g", "h"]

        pd.testing.assert_series_equal(y, y_expect)

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

        pd.testing.assert_frame_equal(
            coef.sort_values(by="feature"), coef_exp.sort_values(by="feature")
        )

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
