import yaml

def load_config(config_file):
    """
    Load configuration from a YAML file.

    Args:
    config_file (str): Path to the YAML configuration file.

    Returns:
    dict: A dictionary containing the configuration parameters.
    """
    try:
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)
        print(f"Configuration loaded successfully from {config_file}")
        return config
    except FileNotFoundError:
        print(f"Error: Configuration file '{config_file}' not found.")
        return None
    except yaml.YAMLError as e:
        print(f"Error parsing YAML file: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred while loading the configuration: {e}")
        return None