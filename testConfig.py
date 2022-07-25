from node import Pipeline
from utils import read_textfile


pipeline = Pipeline('config/config_stage5.yaml', 'test_results/stage4/consensus.txt')
pipeline.run()

# with open('config/template.yaml', 'r') as configfile:
#     config = yaml.load(configfile, Loader=yaml.FullLoader)
#
# items = config['NodeFunction']
# for item in items:
#     print(item)