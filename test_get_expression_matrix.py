from webbrowser import get
from get_expression_matrix import get_expression_matrix
import pandas as pd


def test_alwaysPass():
    assert True


def test_alwaysFail():
    assert False


def test_example_cropped():
    # run script
    get_expression_matrix("examples/cropped")
    test_df = pd.read_excel("examples/cropped/result/test_result.xlsx")
    df = pd.read_excel("examples/cropped/result/result.xlsx")
    assert test_df.equals(df)
