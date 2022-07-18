from node import Pipeline
from utils import read_textfile


pipeline = Pipeline('config/config_stage2.yaml', 'data_test.txt')
data = read_textfile("data_test.txt")
pipeline.run()

# with open('config/template.yaml', 'r') as configfile:
#     config = yaml.load(configfile, Loader=yaml.FullLoader)
#
# items = config['NodeFunction']
# for item in items:
#     print(item)