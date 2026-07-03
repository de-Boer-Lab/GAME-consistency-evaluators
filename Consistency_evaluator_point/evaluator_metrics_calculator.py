'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented
# This evaluator measures reverse-complement (RC) consistency: for each sequence S, 
# the Predictor's prediction for S should track its prediction for RC(S)


import os
import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime, timezone
from scipy.stats import pearsonr

from config import EVALUATOR_NAME

def calculate_pearson_r(predictions_content: str):
    """
    Calculates the Pearson correlation coefficient (r) between sequences and their reverse complements.

    Args:
        predictions_json_path (str): Path to JSON file with predictions
    Returns:
        float : the correlation coefficient r
        0.0   : predictions were returned but one side has zero variance
        None  : disqualified -- no predictions, an error payload, NA/non-numeric
                values. The caller writes this as 'NaN'
    """
    predictions_dict = predictions_content['prediction_tasks'][0]['predictions']

    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        return None
    
    # Flatten {id: value} into a frame, unwrapping single-element lists ([x] -> x).
    rows = []
    for seq_id, value in predictions_dict.items():
        if isinstance(value, list):
            value = value[0] if len(value) > 0 else None
        rows.append((seq_id, value))
    predictions_df = pd.DataFrame(rows, columns=['seq_name', 'predicted_value'])

    # Any NA prediction disqualifies the task (nothing reliable to correlate).
    na_rows = predictions_df[predictions_df['predicted_value'].isna()]
    if not na_rows.empty:
        print("NA values were found in the predictions -> skipping (pearson_r = NaN).")
        print(na_rows)
        return None
    
    # Separate RC and forward predictions based on 'RC' in the sequence name
    # Pairing by key -- not by row position
    is_rc = predictions_df['seq_name'].str.lower().str.endswith('_rc')
    forward_df = predictions_df[~is_rc].copy()
    rc_df = predictions_df[is_rc].copy()
    rc_df['seq_name'] = rc_df['seq_name'].str.slice(stop=-3)  # Remove '_RC' suffix to align with forward_df
    
    merged = forward_df.merge(rc_df, on='seq_name', suffixes=('_fwd', '_rc'))
    
    # Enforce all or nothing. Predictors must score every sequence they were sent
    # or return error -- no cherry-picking allowed. If any sequence is missing a prediction, skip the evaluation.
    # The main evaluator also gates on request-vs-returned counts; this enforces the same policy at the 
    # metric layer and protects against offline re-scoring that bypasses that gate
    fwd_orphans = set(forward_df['seq_name']) - set(merged['seq_name'])
    rc_orphans = set(rc_df['seq_name']) - set(merged['seq_name'])
    if fwd_orphans or rc_orphans:
        print(f"Incomplete prediction set: {len(fwd_orphans)} forward and "
              f"{len(rc_orphans)} RC sequence(s) are missing their partner -> "
              f"refusing to score (pearson_r = NaN)")
        examples = list(fwd_orphans | rc_orphans)[:5]
        if examples:
            print(f"  e.g. missing a partner for: {examples}")
        return None

    if len(merged) < 2:
        print("Fewer than 2 forward/RC pairs -> cannot compute correlation (pearson_r = NaN).")
        return None

    x = pd.to_numeric(merged['predicted_value_fwd'], errors='coerce')
    y = pd.to_numeric(merged['predicted_value_rc'], errors='coerce')
    if x.isna().any() or y.isna().any():
        print("Non-numeric prediction values found -> skipping (pearson_r = NaN).")
        return None

    # Zero variance: predictions exist but one side is constant, so pearsonr is
    # undefined. The framework scores a run that produced predictions as 0.0; NaN is
    # reserved for 'no predictions'.
    if x.std() == 0 or y.std() == 0:
        print("Zero variance in forward or RC predictions -> pearson_r = 0.0")
        return 0.0

    r, _ = pearsonr(x, y)
    pearson_r = 0.0 if np.isnan(r) else float(r)
    print(f"Calculated Pearson r: {pearson_r}")
    return pearson_r


def calculate_and_save_metrics(saved_predictions_path, output_dir):
    """
    Calculates custom evaluation metrics and saves them to CSV files.
    This is the primary function to customize for a new evaluator.
    """
    
    try:
        if not os.path.exists(saved_predictions_path):
            print(f"Predictions file not found: {saved_predictions_path}", file=sys.stderr)
            return

        print("Starting evaluation calculation and saving as CSV.")
        print(f"Using predictions from: {saved_predictions_path}")
        print(f"Summary will be written to: {output_dir}")

        with open(saved_predictions_path, 'r') as f:
            predictions_file_content = json.load(f)

        # Guard the predictor name so a missing key can't AttributeError on .replace().
        predictor_name_received = predictions_file_content.get("predictor_name") or "UnknownPredictor"
        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")

        description = "Consistency point (K562)"
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")

        pearson_r = calculate_pearson_r(predictions_file_content)

        # Task metadata without the bulky predictions payload.
        prediction_task_data_onlyinfo = [
            {k: v for k, v in predictions_file_content["prediction_tasks"][0].items()
             if k != "predictions"}
        ]

        # None (disqualified / no predictions) -> 'NaN'; 0.0 (ran, zero variance) stays 0.0.
        value_str = "NaN" if pearson_r is None else str(pearson_r)

        evaluation_output = pd.DataFrame([{
            'evaluator_name': EVALUATOR_NAME,
            'description': description,
            'predictor_name': predictor_name,
            'time_stamp': timestamp,
            'metric': 'pearson_r',
            'value': value_str,
            'prediction_task(s)_data': prediction_task_data_onlyinfo,
        }])

        summary_filepath = os.path.join(output_dir, f"evaluation_summary_{EVALUATOR_NAME}.csv")
        file_exists = os.path.isfile(summary_filepath)
        evaluation_output.to_csv(summary_filepath, mode='a', sep='\t',
                                 header=(not file_exists), index=False)
        if file_exists:
            print(f"Appended metrics to {summary_filepath}")
        else:
            print(f"Created new metrics file {summary_filepath}")

    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
    
