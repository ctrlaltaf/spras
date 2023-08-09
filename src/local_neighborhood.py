from src.util import prepare_volume, run_container

import pandas as pd

from pathlib import Path

from src.prm import PRM

__all__ = ["LocalNeighborhood"]


class LocalNeighborhood(PRM):
    required_inputs = ["network", "nodes"]

    def generate_inputs(data, filename_map):
        """
        Access fields from the dataset and write the required input files
        @param data: dataset
        @param filename_map: a dict mapping file types in the required_inputs to the filename for that type
        @return:
        """
        for input_type in LocalNeighborhood.required_inputs:
            if input_type not in filename_map:
                raise ValueError(f"{input_type} filename is missing")

        if data.contains_node_columns("prize"):
            # NODEID is always included in the node table
            node_df = data.request_node_columns(["prize"])
        elif data.contains_node_columns(["sources", "targets"]):
            # If there aren't prizes but are sources and targets, make prizes based on them
            node_df = data.request_node_columns(["sources", "targets"])
        else:
            raise ValueError(
                "Omics Integrator 1 requires node prizes or sources and targets"
            )

        print(node_df.to_string())

        node_df.to_csv(
            filename_map["nodes"], index=False, columns=["NODEID"], header=False
        )

        # For now we assume all input networks are undirected until we expand how edge tables work
        edges_df = data.get_interactome()
        print(edges_df)
        edges_df.to_csv(
            filename_map["network"],
            sep="|",
            index=False,
            columns=["Interactor1", "Interactor2"],
            header=False,
        )

    @staticmethod
    def run(nodes=None, network=None, output_file=None, singularity=False):
        """
        Run PathLinker with Docker
        @param nodes:  input node types with sources and targets (required)
        @param network:  input network file (required)
        @param output_file: path to the output pathway file (required)
        @param singularity: if True, run using the Singularity container instead of the Docker container
        """
        # Add additional parameter validation
        # Do not require k
        # Use the PathLinker default
        # Could consider setting the default here instead
        if not nodes or not network or not output_file:
            raise ValueError("Required LocalNeighborhood arguments are missing")

        work_dir = "/spras"

        # Each volume is a tuple (src, dest)
        volumes = list()

        bind_path, node_file = prepare_volume(nodes, work_dir)
        volumes.append(bind_path)

        bind_path, network_file = prepare_volume(network, work_dir)
        volumes.append(bind_path)

        bind_path, bound_output_file = prepare_volume(output_file, work_dir)
        volumes.append(bind_path)

        command = [
            "python",
            "/LocalNeighborhood/local_neighborhood.py",
            "--network",
            network_file,
            "--nodes",
            node_file,
            "--output",
            bound_output_file,
        ]

        print(
            "Running Local Neighborhood with arguments: {}".format(" ".join(command)),
            flush=True,
        )

        # TODO consider making this a string in the config file instead of a Boolean
        container_framework = "singularity" if singularity else "docker"
        out = run_container(
            container_framework,
            "ctrlaltaf/local-neighborhood",
            command,
            volumes,
            work_dir,
        )
        print(out)

    @staticmethod
    def parse_output(raw_pathway_file, standardized_pathway_file):
        """
        Convert a predicted pathway into the universal format
        @param raw_pathway_file: pathway file produced by an algorithm's run function
        @param standardized_pathway_file: the same pathway written in the universal format
        """
        # Questions: should there be a header/optional columns?
        # What about multiple raw_pathway_files
        df = pd.read_csv(raw_pathway_file, sep="|")
        df.insert(2, "rank", 1)  # Add a constant rank of 1
        df.to_csv(standardized_pathway_file, sep="\t", header=False, index=False)
