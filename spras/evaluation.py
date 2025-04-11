import os
import pickle as pkl
from pathlib import Path
from typing import Dict, Iterable
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.metrics import precision_score, recall_score


class Evaluation:
    NODE_ID = "NODEID"

    def __init__(self, gold_standard_dict: Dict):
        self.label = None
        self.datasets = None
        self.node_table = None
        self.edge_table = None # TODO: later iteration
        self.load_files_from_dict(gold_standard_dict)
        return

    @staticmethod
    def merge_gold_standard_input(gs_dict, gs_file):
        """
        Merge files listed for this gold standard dataset and write the dataset to disk
        @param gs_dict: gold standard dataset to process
        @param gs_file: output filename
        """
        gs_dataset = Evaluation(gs_dict)
        gs_dataset.to_file(gs_file)

    def to_file(self, file_name):
        """
        Saves gold standard object to pickle file
        """
        with open(file_name, "wb") as f:
            pkl.dump(self, f)

    @staticmethod
    def from_file(file_name):
        """
        Loads gold standard object from a pickle file.
        Usage: gold_standard = Evaluation.from_file(pickle_file)
        """
        with open(file_name, "rb") as f:
            return pkl.load(f)

    def load_files_from_dict(self, gold_standard_dict: Dict):
        """
        Loads gold standard files from gold_standard_dict, which is one gold standard dataset
        dictionary from the list in the config file with the fields in the config file.
        Populates node_table.

        node_table is a single column of nodes pandas table.

        returns: none
        """
        self.label = gold_standard_dict["label"]  # cannot be empty, will break with a NoneType exception
        self.datasets = gold_standard_dict["dataset_labels"]  # can be empty, snakemake will not run evaluation due to dataset_gold_standard_pairs in snakemake file

        # cannot be empty, snakemake will run evaluation even if empty
        node_data_files = gold_standard_dict["node_files"][0]  # TODO: single file for now

        data_loc = gold_standard_dict["data_dir"]

        single_node_table = pd.read_table(os.path.join(data_loc, node_data_files), header=None)
        single_node_table.columns = [self.NODE_ID]
        self.node_table = single_node_table

        # TODO: are we allowing multiple node files or single in node_files for gs
        # if yes, a for loop is needed

        # TODO: later iteration - chose between node and edge file, or allow both
        edge_data_files = gold_standard_dict["edge_files"][0]
        single_edge_table = pd.read_table(os.path.join(data_loc, edge_data_files), header=None)
        self.edge_table = single_edge_table

    @staticmethod
    def precision(file_paths: Iterable[Path], node_table: pd.DataFrame, output_file: str):
        """
        Takes in file paths for a specific dataset and an associated gold standard node table.
        Calculates precision for each pathway file
        Returns output back to output_file
        @param file_paths: file paths of pathway reconstruction algorithm outputs
        @param node_table: the gold standard nodes
        @param output_file: the filename to save the precision of each pathway
        """
        y_true = set(node_table['NODEID'])
        results = []

        for file in file_paths:
            df = pd.read_table(file, sep="\t", header=0, usecols=["Node1", "Node2"])
            y_pred = set(df['Node1']).union(set(df['Node2']))
            all_nodes = y_true.union(y_pred)
            y_true_binary = [1 if node in y_true else 0 for node in all_nodes]
            y_pred_binary = [1 if node in y_pred else 0 for node in all_nodes]

            # default to 0.0 if there is a divide by 0 error
            precision = precision_score(y_true_binary, y_pred_binary, zero_division=0.0)

            results.append({"Pathway": file, "Precision": precision})

        precision_df = pd.DataFrame(results)
        precision_df.to_csv(output_file, sep="\t", index=False)

    # will probably need to split the precision and recall
    @staticmethod
    def precision_and_recall_edge(file_paths: Iterable[Path], edge_table: pd.DataFrame, algorithms: list, output_file: str, output_png:str=None):
        """
        Takes in file paths for a specific dataset and an associated gold standard edge table.
        Calculates precision and recall for each pathway file
        Returns output back to output_file
        @param file_paths: file paths of pathway reconstruction algorithm outputs
        @param edge_table: the gold standard edges
        @param algorithms: list of algorithms used in current run of SPRAS
        @param output_file: the filename to save the precision and recall of each pathway
        @param output_png (optional): the filename to plot the precision and recall of each pathway (not a PRC)
        """
        print("EDGE PR")
        algorithm_directionality_dict = {"pathlinker": "U", "omicsintegrator1": "D"}

        y_true = set()
        for row in edge_table.itertuples():
            y_true.add((row[1], row[2]))
        results = []
        for file in file_paths:
            algorithm = file.split("/")[1].split("-")[1]
            df = pd.read_table(file, sep="\t", header=0, usecols=["Node1", "Node2"])
            y_pred = set()
            for row in df.itertuples():
                y_pred.add((row[1], row[2]))
            all_edges = set(y_true.union(y_pred))

            y_true_binary = []
            y_pred_binary = []
            if algorithm_directionality_dict[algorithm] == "U":
                y_true_binary = [1 if (edge[0], edge[1]) in y_true or (edge[1], edge[0]) in y_true else 0 for edge in all_edges]
                y_pred_binary = [1 if (edge[0], edge[1]) in y_pred or (edge[1], edge[0]) in y_pred else 0 for edge in all_edges]
            else:
                y_true_binary = [1 if (edge[0], edge[1]) in y_true else 0 for edge in all_edges]
                y_pred_binary = [1 if (edge[0], edge[1]) in y_pred else 0 for edge in all_edges]

            precision = precision_score(y_true_binary, y_pred_binary, zero_division=0.0)
            recall = recall_score(y_true_binary, y_pred_binary, zero_division=0.0)
            results.append({"Pathway": file, "Precision": precision, "Recall": recall})

        pr_df = pd.DataFrame(results)
        pr_df.sort_values(by=["Recall", "Pathway"], axis=0, ascending=True, inplace=True)
        pr_df.to_csv(output_file, sep="\t", index=False)

        num_of_algorithms_used = 0
        if output_png is not None:
            if not pr_df.empty:
                plt.figure(figsize=(8, 6))
                # plot a line per algorithm
                for algorithm in algorithms: #TODO I think there is a better way than doing this; using split on the filepaths doesn't work bc it is not adaptable
                    subset = pr_df[pr_df["Pathway"].str.contains(algorithm)]
                    if not subset.empty:
                        plt.plot(
                            subset["Recall"],
                            subset["Precision"],
                            marker='o',
                            linestyle='-',
                            label=f"{algorithm}"
                        )
                        num_of_algorithms_used += 1

                # plot overall precision and recall from all the algorithms
                if num_of_algorithms_used > 1:
                    plt.plot(pr_df["Recall"], pr_df["Precision"], marker='o', linestyle='-', color='b', label="Overall Precision-Recall")

                plt.xlabel("Recall")
                plt.ylabel("Precision")
                plt.title(f"Precision and Recall Plot")
                plt.legend()
                plt.grid(True)
                plt.savefig(output_png)
            else:
                plt.figure()
                plt.plot([], [])
                plt.title("Empty Pathway Files")
                plt.savefig(output_png)


    @staticmethod
    def jaccard_edge_heatmap(file_paths: Iterable[Path], edge_table: pd.DataFrame, output_png:str=None):
        """
        Takes in file paths for a specific dataset and an associated gold standard edge table.
        Generates a jaccard index heatmap image that compares all the edge similarity between each dataset and the gold standard
        Returns output back to output_png
        @param file_paths: file paths of pathway reconstruction algorithm outputs
        @param edge_table: the gold standard edges
        @param output_png (optional): the filename to plot the heatmap (not a PRC)
        """
        print("JACCARD INDEX")
        algorithm_directionality_dict = {"pathlinker": "U", "omicsintegrator1": "D"}

        gs_edges = set()
        for row in edge_table.itertuples():
            gs_edges.add((row[1], row[2]))
        
        # calculate all the jaccard edge index for each method against the gold standard
        jaccard_edge_indices_list = []
        algorithms = []
        for file in file_paths:
            algorithm = file.split("/")[1].split("-")[1]
            df = pd.read_table(file, sep="\t", header=0, usecols=["Node1", "Node2"])
            method_edges = set()
            for row in df.itertuples():
                if algorithm_directionality_dict[algorithm] == "U":
                    method_edges.add((row[1], row[2]))
                    method_edges.add((row[2], row[1]))
                else:
                    method_edges.add((row[1], row[2]))
            edge_union = gs_edges | method_edges
            edge_intersection = gs_edges & method_edges
            jaccard_edge_index = len(edge_intersection) / len(edge_union)
            jaccard_edge_indices_list.append(float(jaccard_edge_index))
            algorithms.append(f"{file.split('/')[1].split('-')[1]}-{file.split('/')[1].split('-')[3]}")

        jaccard_edge_indices = np.asanyarray([jaccard_edge_indices_list])

        plt.figure(figsize=(10, 8))
        sns.heatmap(
            jaccard_edge_indices,
            annot=True,
            cmap="viridis",
            xticklabels=algorithms,
            yticklabels=[""],
        )
        plt.xlabel("Algorithms")
        plt.ylabel("Pathways")
        plt.title("Jaccard Index Edge Heatmap")
        plt.tick_params(axis='x', which='major', labelsize=7.5)
        plt.savefig(output_png, format="png", dpi=300)