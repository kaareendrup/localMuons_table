import json

def format_json_indents(input_path):
    # Load the JSON file
    with open(input_path, "r") as f:
        data = json.load(f)

    output_path = input_path.replace(".json", "_formatted.json")

    # Save the formatted JSON with indents
    with open(output_path, "w") as f:
        json.dump(data, f, indent=4)

    print(f"Indented JSON saved to {output_path}")

input_file = "config/configuration.json"
format_json_indents(input_file)