from node import Pipeline
from utils import read_textfile


pipeline = Pipeline('config/config_stage2.yaml', 'test_results/stage2/output.fastq')
pipeline.run()

# with open('config/template.yaml', 'r') as configfile:
#     config = yaml.load(configfile, Loader=yaml.FullLoader)
#
# items = config['NodeFunction']
# for item in items:
#     print(item)