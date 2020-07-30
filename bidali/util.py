# -*- coding: utf-8 -*-
"""bidali utilities

Functions for transforming pd.DataFrame and other generic data processing steps.id
"""
import pandas as pd

def unfoldDFlistColumn(df,listColName):
    """Unfold a pd.DataFrame

    A dataframe containing a column with lists as values will be unfolded to
    a dataframe with a new row for every value in the list.

    Args:
        df (pd.DataFrame): DataFrame to unfold.
        listColName (str): Name of the column that has list values.
    """
    unfoldedListCol = df.apply(
        lambda x: pd.Series(x[listColName]),axis=1
    ).stack().reset_index(level=1, drop=True)
    unfoldedListCol.name = listColName
    unfoldedDF = df.drop(listColName, axis=1).join(unfoldedListCol)
    return unfoldedDF
