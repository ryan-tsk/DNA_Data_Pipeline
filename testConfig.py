from node import Pipeline

pipeline = Pipeline('archived/config/config_stage4.yaml', 'test_data/short_data.txt')
pipeline.run()
