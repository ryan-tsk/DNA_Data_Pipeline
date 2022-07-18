from node import Pipeline
from utils import read_textfile


pipeline = Pipeline('config/config_stage1.yaml')
data = read_textfile("data_test.txt")
pipeline.run(data)

# with open('config/template.yaml', 'r') as configfile:
#     config = yaml.load(configfile, Loader=yaml.FullLoader)
#
# items = config['NodeFunction']
# for item in items:
#     print(item)