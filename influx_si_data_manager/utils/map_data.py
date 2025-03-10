import pandas as pd

CONNECTOR_NAME = "influx_si_data_connector"

def map_data(mapping_file, data, from_tool):
    """
    Replace cell values in dataset by following indications in mapping file
    :param mapping_file: Path to the mapping file
    :param data: dataset to modify
    :param from_tool: tool name to filter mapping entries
    :return: Mapped data
    """
    # Read the mapping file and filter for relevant entries
    mapping = pd.read_csv(mapping_file, sep="\t")
    filtered_mapping = mapping[(mapping["Connector"] == CONNECTOR_NAME) & 
                              (mapping["From_tool"] == from_tool)]
    
    # Replace values in each column
    for column in filtered_mapping["Column"].unique():
        replacements = filtered_mapping[filtered_mapping["Column"] == column]
        replace_dict = dict(zip(replacements["Replace"], replacements["By"]))
        data[column] = data[column].replace(replace_dict)
    
    return data
            
            

if __name__ == "__main__":
    
    isocor_data = pd.read_csv("../test_data/isocor_results.tabular", sep="\t")
    print(isocor_data)
    datasets = {"isocor": isocor_data}
    replaced = map_data(r"../test_data/mapping.txt", data=isocor_data, from_tool="isocor")
    print(replaced)