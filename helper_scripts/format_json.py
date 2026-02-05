import json

def json_to_string(input_path, key_name="embedded"):
    # Load the inner JSON file
    with open(input_path, "r") as f:
        inner_json = json.load(f)

    output_path = input_path.replace(".json", "_formatted.json")

    # Convert to a compact escaped string
    json_string = json.dumps(inner_json)

    # Build the outer JSON
    outer = {key_name: json_string}

    # Save the result
    with open(output_path, "w") as f:
        json.dump(outer, f, indent=4)

    print(f"Formatted JSON saved to {output_path}")

path = "cuts/"
# name = "custom_cuts_global_split"
# name = "custom_cuts_global"
name = "custom_cuts_standalone"

json_to_string(f"{path}{name}.json", key_name="cfgMuonCutsJSON")
