"""
Methods and functions to facilitate machine learning analyses of cell health data
"""

import os
import numpy as np
import pandas as pd

from sklearn.cluster import KMeans
from sklearn.preprocessing import minmax_scale
from sklearn.metrics import (
    mean_squared_error,
    r2_score,
    roc_auc_score,
    roc_curve,
    precision_recall_curve,
    average_precision_score,
)
from joblib import dump, load
from dask_ml.model_selection import GridSearchCV


class CellHealthPredict:
    def __init__(
        self,
        x_df,
        y_df,
        parameters,
        estimator,
        n_folds=5,
        cv_scoring="neg_mean_squared_error",
        return_train_score=True,
        shuffle=False,
    ):
        if shuffle:
            self.x_df = x_df.apply(shuffle_columns, axis=1, result_type="broadcast")
            self.shuffle_key = "shuffle_true"
        else:
            self.x_df = x_df
            self.shuffle_key = "shuffle_false"
        self.y_df = y_df
        self.parameters = parameters
        self.estimator = estimator
        self.named_step = self.estimator.steps[0][0]
        self.n_folds = n_folds
        self.cv_scoring = cv_scoring
        self.return_train_score = return_train_score
        self.profile_ids = self.x_df.index.tolist()
        self.is_fit = False
        self.kmeans_fit = False
        self.cv_pipeline = GridSearchCV(
            estimator=self.estimator,
            param_grid=self.parameters,
            n_jobs=-1,
            cv=self.n_folds,
            scoring=self.cv_scoring,
            return_train_score=self.return_train_score,
        )

    def fit_cell_health_target(self, target, y_transform="raw", binarize_fit="kmeans"):
        """
        Fit the cross validation pipeline to obtain an optimal model

        Arguments:
        target - a string that maps to a specific column in the y dataframe
        y_transform - a string indicating how to transform the target variable
                      [default: "raw"; can be "raw", "zero-one", or "binarize"]

        Output:
        Model is fit internally and can be accessed through self
        """
        # Assert that the feature exists in the initialized target DataFrame
        self.target = target
        assert self.target in self.y_df.columns.tolist(), "Feature must exist in y_df"

        # Subset y_df to the given target
        self.y = self.y_df.loc[:, self.target]

        # Assert that the user specified y_transform is supported
        one_of_y_transform = ["raw", "zero-one", "binarize"]
        assert (
            y_transform in one_of_y_transform
        ), "'{}' is not supported try one of: {}".format(
            y_transform, ", ".join(one_of_y_transform)
        )
        self.y_transform = y_transform

        # Realign the input matrices
        self.x_realigned_df, self.y_scaled_df, self.n_samples_removed = self.realign_missing_data(
            self.x_df, self.y
        )

        # Update the profile IDs
        self.profile_ids = self.x_realigned_df.index

        # Transform the target variable
        self.y_scaled_df = self.recode_y(y=self.y_scaled_df, binarize_fit=binarize_fit)

        # The model can proceed only if there are enough representation across classes
        if y_transform == "binarize":
            class_count = self.y_scaled_df.value_counts()
            if (class_count > self.n_folds).sum() < 2:
                return False
            self.min_class_count = class_count.min()

        # Fit the cross validation pipeline object
        self.pipeline_fit = self.cv_pipeline.fit(
            X=self.x_realigned_df, y=self.y_scaled_df
        )
        self.is_fit = True
        return True

    def realign_missing_data(self, x_df, y):
        """
        Given two input dataframes, remove missing values and realign

        Arguments:
        x_df - input DataFrame of features
        y - input Series of target variables

        Output:
        The realigned X and y matrices and the number of samples removed
        """
        # Sometimes the target includes missing values, make sure to check and remove
        pd.testing.assert_index_equal(x_df.index, y.index, check_names=False)

        # How many samples are removed
        num_samples_removed = y.isna().sum()

        # Drop these samples and realign to x matrix
        y = y.dropna()
        x_df = x_df.loc[y.index, :]

        # Assert that the matrices have equivalent indices
        pd.testing.assert_index_equal(x_df.index, y.index, check_names=False)

        return x_df, y, num_samples_removed

    def recode_y(self, y, target="None", y_transform="None", binarize_fit="kmeans"):
        """
        Recode y using a specific transformation. Note: if no arguments are passed, the
        function will use elements stored in self.

        Arguments:
        y - pandas Series indicating the status of the target variable
        target - "None" or string indicating target variable in y [Default: "None"]
        y_transform - "None" or string indicating transform operation [Default: "None"]

        Output:
        The transformed y variables
        """
        if target == "None":
            target = self.target

        if y_transform == "None":
            y_transform = self.y_transform

        if isinstance(y, pd.DataFrame):
            y = y.loc[:, target]

        if y_transform == "zero-one":
            y_scale = minmax_scale(y)
        elif y_transform == "z-score":
            y_scale = scale(y)
        elif y_transform == "binarize":
            if binarize_fit == "kmeans":
                if not self.kmeans_fit:
                    self.kmeans = KMeans(n_clusters=2, random_state=0).fit(
                        np.reshape(y.values, (-1, 1))
                    )
                    self.kmeans_fit = True
                    y_scale = self.kmeans.labels_
                else:
                    y_scale = self.kmeans.predict(np.reshape(y.values, (-1, 1)))
            else:
                feature_median = y.median()
                y_neg = y < feature_median
                y_pos = y >= feature_median
                y_scale = y.copy()

                y_scale.loc[y_neg] = 0
                y_scale.loc[y_pos] = 1
                y_scale = y_scale.astype(int)
        else:
            y_scale = y

        y_scale = pd.Series(y_scale, name=target)
        y_scale.index = y.index

        return y_scale

    def set_y_transform(self, y_transform):
        self.y_transform = y_transform

    def set_target(self, target):
        self.target = target

    def get_cv_results(self):
        assert (
            self.is_fit
        ), "The model is not yet fit! Run fit_cell_health_target() first"

        cv_results_df = pd.concat(
            [
                pd.DataFrame(self.pipeline_fit.cv_results_).drop(
                    "params", axis="columns"
                ),
                pd.DataFrame.from_records(self.pipeline_fit.cv_results_["params"]),
            ],
            axis="columns",
        )

        cv_results_df = (
            pd.pivot_table(
                cv_results_df,
                values="mean_test_score",
                index="{}__l1_ratio".format(self.named_step),
                columns="{}__alpha".format(self.named_step),
            )
            .stack()
            .reset_index()
            .assign(
                target=self.target,
                shuffle=self.shuffle_key,
                y_transform=self.y_transform,
            )
        )

        cv_results_df.columns = [
            "l1_ratio",
            "alpha",
            "mean_test_score",
            "target",
            "shuffle",
            "y_transform",
        ]

        return cv_results_df

    def predict(self, x, decision_function=False):
        """
        If decision_function, output continuous label
        """
        assert (
            self.is_fit
        ), "The model is not yet fit! Run fit_cell_health_target() first"
        if decision_function:
            assert (
                self.y_transform == "binarize"
            ), "decision_function only for classifiers not regressors"
            y_pred = self.pipeline_fit.decision_function(x)
        else:
            y_pred = self.pipeline_fit.predict(x)

        return pd.Series(y_pred, name=self.target)

    def get_coefficients(self, save_model=False, model_file="tmp.joblib"):
        """
        Extract coefficients from the model

        Arguments:
        save_model - boolean if the model should be saved or not
        model_file - string indicating where the model should be saved

        Output: Coefficients of the model
        """
        # Save classifier coefficients
        final_pipeline = self.cv_pipeline.best_estimator_
        final_classifier = final_pipeline.named_steps[self.named_step]

        if self.named_step == "classify":
            weights = final_classifier.coef_[0]
        else:
            weights = final_classifier.coef_

        coef_df = pd.DataFrame.from_dict(
            {"feature": self.x_df.columns, "weight": weights}
        )

        coef_df = coef_df.assign(
            abs_weight=coef_df.weight.abs(),
            target=self.target,
            y_transform=self.y_transform,
            shuffle=self.shuffle_key,
        )
        coef_df = coef_df.sort_values("abs_weight", ascending=False).reset_index(
            drop=True
        )

        if save_model:
            dump(final_classifier, model_file)

        return coef_df

    def get_performance(
        self,
        decision_function=False,
        x_test=None,
        y_test=None,
        return_y=False,
        binarize_fit="kmeans",
    ):
        """
        Get the classifier or regression performance of the fit classifier

        Arguments:
        decision_function - boolean if the output should return continuous output [For "binarize" only]
        x_test - "None" or pandas DataFrame. If None, get training performance [Default: "None"]
        y_test - "None" or pandas DataFrame indicating status labels. [Default: "None"]
        return_y - boolean if the transformed y variable should be output [Default: False]
        binarize_fit - string indicating how the binary recoding is performed

        Output:
        Performanc metrics for regression or classification (depends on y_transform)
        """

        # Make sure the model has been fit first!
        if not self.is_fit:
            ValueError("The model is not yet fit! Run fit_cell_health_target() first")

        # If the input is a dataframe, make sure columns are aligned to training data
        if isinstance(x_test, pd.DataFrame):
            assert (
                x_test.columns.tolist() == self.x_realigned_df.columns.tolist()
            ), "Features are not aligned!"

            pd.testing.assert_index_equal(x_test.index, y_test.index, check_names=False)

            # Make sure missing data is removed
            y_true = y_test.loc[:, self.target]
            x_test, y_true, self.n_test_samples_removed = self.realign_missing_data(
                x_test, y_true
            )

            # Transform the target variable
            y_true = self.recode_y(y=y_true, binarize_fit=binarize_fit)

            # Get predicted values from trained model
            y_pred = self.predict(x_test, decision_function=decision_function)
            data_fit_type = "test"
            profile_ids = x_test.index.tolist()

        else:
            y_pred = self.predict(
                self.x_realigned_df, decision_function=decision_function
            )
            y_true = self.y_scaled_df
            data_fit_type = "train"
            profile_ids = self.x_realigned_df.index.tolist()

        if self.y_transform != "binarize":
            mse = mean_squared_error(y_true, y_pred)
            r_two = r2_score(y_true, y_pred)

            mse_df = pd.DataFrame([mse], columns=["mse"])
            mse_df = mse_df.assign(
                metric="mse",
                target=self.target,
                data_fit=data_fit_type,
                shuffle=self.shuffle_key,
                y_transform=self.y_transform,
            )

            r_two_df = pd.DataFrame([r_two], columns=["mse"])
            r_two_df = r_two_df.assign(
                metric="r_two",
                target=self.target,
                data_fit=data_fit_type,
                shuffle=self.shuffle_key,
                y_transform=self.y_transform,
            )

            output = [mse_df, r_two_df]

        else:
            auroc_weighted = roc_auc_score(y_true, y_pred, average="weighted")
            aupr_weighted = average_precision_score(y_true, y_pred, average="weighted")

            roc_columns = ["fpr", "tpr", "threshold"]
            pr_columns = ["precision", "recall", "threshold"]

            roc_items = zip(roc_columns, roc_curve(y_true, y_pred))
            roc_df = pd.DataFrame.from_dict(dict(roc_items))
            roc_df = roc_df.assign(
                metric="roc",
                target=self.target,
                auc=auroc_weighted,
                data_fit=data_fit_type,
                shuffle=self.shuffle_key,
                y_transform=self.y_transform,
                min_class_count=self.min_class_count,
            )

            prec, rec, thresh = precision_recall_curve(y_true, y_pred)
            pr_df = pd.DataFrame.from_records([prec, rec]).transpose()
            pr_df = pd.concat([pr_df, pd.Series(thresh)], ignore_index=True, axis=1)
            pr_df.columns = pr_columns
            pr_df = pr_df.assign(
                metric="aupr",
                target=self.target,
                auc=aupr_weighted,
                data_fit=data_fit_type,
                shuffle=self.shuffle_key,
                y_transform=self.y_transform,
                min_class_count=self.min_class_count,
            )

            output = [roc_df, pr_df]

        if return_y:

            y_true = (
                pd.DataFrame(y_true)
                .rename({self.target: "recode_target_value"}, axis="columns")
                .assign(
                    target=self.target,
                    data_type=data_fit_type,
                    shuffle=self.shuffle_key,
                    y_transform=self.y_transform,
                    y_type="y_true",
                )
            )
            y_true.index.name = "Metadata_profile_id"
            y_true = y_true.reset_index()

            y_pred = (
                pd.DataFrame(y_pred)
                .rename({self.target: "recode_target_value"}, axis="columns")
                .assign(
                    target=self.target,
                    data_type=data_fit_type,
                    shuffle=self.shuffle_key,
                    y_transform=self.y_transform,
                    y_type="y_pred",
                )
            )
            y_pred.index = profile_ids
            y_pred.index.name = "Metadata_profile_id"
            y_pred = y_pred.reset_index()

            output += [y_true, y_pred]

        return output


def load_train_test(data_dir="data", drop_metadata=False, output_metadata_only=False):
    """
    Load training and testing data

    Arguments:
    data_dir - a string indicating the location of the files
    drop_metadata - boolean if the metadata should be striped before input
                    [default: False]
    output_metadata_only - boolean if metadata columns are only output
                           [default: False]

    Output:
    A list storing training and testing x and y matrices
    """

    x_train_file = os.path.join(data_dir, "x_train.tsv.gz")
    y_train_file = os.path.join(data_dir, "y_train.tsv.gz")
    x_test_file = os.path.join(data_dir, "x_test.tsv.gz")
    y_test_file = os.path.join(data_dir, "y_test.tsv.gz")

    x_train_df = pd.read_csv(x_train_file, index_col=0, sep="\t")
    y_train_df = pd.read_csv(y_train_file, index_col=0, sep="\t")
    x_test_df = pd.read_csv(x_test_file, index_col=0, sep="\t")
    y_test_df = pd.read_csv(y_test_file, index_col=0, sep="\t")

    x_train_metadata = x_train_df.columns.str.startswith("Metadata_")
    x_test_metadata = x_test_df.columns.str.startswith("Metadata_")

    y_train_metadata = y_train_df.columns.str.startswith("Metadata_")
    y_test_metadata = y_test_df.columns.str.startswith("Metadata_")

    if drop_metadata:
        x_train_df = x_train_df.loc[:, ~x_train_metadata]
        x_test_df = x_test_df.loc[:, ~x_test_metadata]

        y_train_df = y_train_df.loc[:, ~y_train_metadata]
        y_test_df = y_test_df.loc[:, ~y_test_metadata]

    elif output_metadata_only:
        x_train_df = x_train_df.loc[:, x_train_metadata]
        x_test_df = x_test_df.loc[:, x_test_metadata]

        y_train_df = y_train_df.loc[:, y_train_metadata]
        y_test_df = y_test_df.loc[:, y_test_metadata]

    return x_train_df, x_test_df, y_train_df, y_test_df


def shuffle_columns(cp_feature):
    """
    To be used in an `apply` pandas func to shuffle columns around a datafame
    Import only
    """
    import numpy as np

    return np.random.permutation(cp_feature.tolist())
