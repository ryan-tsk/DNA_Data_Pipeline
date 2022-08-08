import argparse

from node import Pipeline

parser = argparse.ArgumentParser(description='DNA Data Pipeline')
parser.add_argument('-i', '--input', help='Input file to be processed. Optional if specified in config', default=None)
parser.add_argument('-s', '--stages', help='Set to true if intermediate results wanted. Optional.', default=False)
parser.add_argument('-c', '--config', required=True, help='Configuration file for pipeline. Required.')
parser.add_argument('-o', '--output',
                    help='Output file for pipeline. Defaults to results.txt in current folder. Optional',
                    default='results.txt')

args = parser.parse_args()
pipeline = Pipeline(config_path=args.config, input_path=args.input, output_path=args.output, stages=args.stages)
