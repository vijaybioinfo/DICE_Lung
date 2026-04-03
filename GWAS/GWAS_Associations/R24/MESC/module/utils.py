import os
import pandas as pd


def read_file(file: str):
    return pd.read_table(file)

def get_filename(path: str) -> str:
    return os.path.basename(path)

def trim_file_extension(filename: str):
    if filename.endswith(".all.h2med"):
        return filename.replace(".all.h2med", "")

    if filename.endswith(".categories.h2med"):
        return filename.replace(".categories.h2med", "")

    return os.path.splitext(filename)[0]

def mkdir(path: str) -> None:
    os.makedirs(path, exist_ok=True)

def parent_directory(path: str) -> str:
    return os.path.abspath(os.path.join(path, os.pardir))

def merge_files(files: list[str], colname: str):
    combined = []
    for file in files:
        data = read_file(file)
        data[colname] = trim_file_extension(get_filename(file))
        combined.append(data)

    combined = pd.concat(combined).reset_index(drop=True)
    return combined