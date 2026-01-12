'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
import itertools
from datetime import datetime, timezone
from scipy.stats import pearsonr

from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH, output_filename_base

def calculate_pearson_r(predictions_content: str):
    """
    Calculates the Pearson correlation coefficient (r) between sequences and their reverse complements.

    Args:
        predictions_json_path (str): Path to JSON file with predictions
    Returns:
        float: The Pearson correlation coefficient (r), or None if calculation isn't possible.
    """
    predictions_dict = predictions_content['prediction_tasks'][0]['predictions']

    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        return None

    # Create DataFrame from Predictions
    predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=['id_column', 'Predicted_Value'])
    predictions_df['Predicted_Value'] = predictions_df['Predicted_Value'].apply(
        lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
        )
    
    #check here is there is NA is any of the prediction values
    na_rows = predictions_df[predictions_df['Predicted_Value'].isna()]
    if not na_rows.empty:
        print("NA values were found in the predictions, skipping evaluation")
        print(na_rows)
        return None

    RC_predictions = predictions_df[predictions_df['id_column'].str.contains('RC', case=False)]
    forward_predictions = predictions_df[~predictions_df['id_column'].str.contains('RC', case=False)]

    print(RC_predictions)
    print(forward_predictions)
    r, _ = pearsonr(forward_predictions['Predicted_Value'], RC_predictions['Predicted_Value'])
    print(f"Calculated Pearson r: {r}") 
    return r


def calculate_and_save_metrics(saved_predictions_path, output_dir):
    """
    Calculates custom evaluation metrics and saves them to CSV files.
    This is the primary function to customize for a new evaluator.
    """
    print("----- Starting Fake Evaluation Calculation and Saving as CSV -----")
    
    try:
        # Correlation calculation
        # NOTE: Every evaluator will do this slightly differently depending on how the data is presented  
        if os.path.exists(saved_predictions_path):
            print("----- Starting Evaluation Calculation and Saving as CSV -----")

            print(f"Using predictions from: {saved_predictions_path}")
            print(f"Correlation metadata will be saved in {output_dir}")
        
        # Now load predictions
        with open(saved_predictions_path, 'r') as f:
            predictions_file_content = json.load(f)

        # Extract Predictor Name
        predictor_name_base = predictions_file_content.get("predictor_name", None) # Resort to None if predictor name is not available

        # ADDITION: Construct file name after receiving predictor_name
        predictor_name_received = predictions_file_content.get("predictor_name", None)
        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")
        output_filename = f"{output_filename_base}_from_{predictor_name}.csv"
  
        description = "Consistency point (K562)"
        # Get UTC timestamp for predictor_nam
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
        # Compute the full RETURN_FILE_PATH using the provided output directory

        RETURN_FILE_PATH = os.path.join(output_dir, output_filename)
        print(f"Will save predictions to: {RETURN_FILE_PATH}")
        
        pearson_r = calculate_pearson_r(predictions_file_content)
        prediction_task_data_onlyinfo = [{k: v for k, v in predictions_file_content["prediction_tasks"][0].items() if k != "predictions"}]

        #add code to create the output file
        evaluation_output = pd.DataFrame([{'Evaluator_name': EVALUATOR_NAME, 'Description': description, 'Predictor_name': predictor_name,  'Time_stamp': timestamp, 'Metric': 'pearson_r', 'Value': str(pearson_r), 'Prediction_task(s)_data': prediction_task_data_onlyinfo}])
        evaluation_output.to_csv(RETURN_FILE_PATH , sep = "\t")

    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
