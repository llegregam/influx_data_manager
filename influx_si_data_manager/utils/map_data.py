import pandas as pd

CONNECTOR_NAME = "influx_si_data_connector"

def map_data(mapping_file, data):
    """
    Replace cell values in dataset by following indications in mapping file
    :param mapping_file: Path to the mapping file
    :param dataset: dataset to modify
    :return: Mapped data
    """

    # Read the mapping file
    mapping = pd.read_csv(mapping_file, sep="\t")
    mapping = mapping.loc[mapping["Connector"] == CONNECTOR_NAME, :]

    for key in data:
        current_tool_view = mapping.loc[mapping["From_tool"] == key]
        for column in current_tool_view["Column"].unique():
            current_column_view = current_tool_view.loc[current_tool_view["Column"] == column]
            to_replace = current_column_view[["Replace", "By"]].set_index("Replace").to_dict()["By"]
            data[key][column] = data[key][column].replace(to_replace)
    
    return data
            
            

if __name__ == "__main__":
    
    isocor_data = pd.read_csv("../test_data/isocor_results.tabular", sep="\t")
    datasets = {"isocor": isocor_data}
    replaced = map_data(r"../test_data/mapping.txt", data=datasets)
