'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator does this slightly differently depending on how the data is
# presented. This evaluator measures reverse-complement (RC) consistency for TRACK
# predictions: for each sequence S, the Predictor returns a per-position track, and the
# track for S should match the track for RC(S) once the RC track is reversed to put it
# back into S's coordinate frame. We report the mean per-sequence Pearson r

import os
import sys
import json
import pandas as pd
import numpy as np
from datetime import datetime, timezone
from scipy.stats import pearsonr

from config import EVALUATOR_NAME


def calculate_pearson_r_track(predictions_content: str):
    """
    Mean per-sequence Pearson r between each forward track and its reverse-complement
    partner's track (the RC track is reversed to align position-by-position).

    Returns:
        float : mean the correlation coefficient r, across all forward/RC pairs
        0.0   : predictions were returned but one side has zero variance
        None  : disqualified -- no predictions, an error payload, NA/non-numeric
                values. The caller writes this as 'NaN'
    """
    predictions_dict = predictions_content['prediction_tasks'][0]['predictions']
    
    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        return None
    
    # Create DataFrame from Predictions
    predictions_df = pd.DataFrame(
        list(predictions_dict.items()), columns=['seq_name', 'predicted_value']
    )
    
    # Catches a whole prediction being null. NaNs inside a track are checked per-pair
    # below, since a scalar .isna() does not look inside a list
    if predictions_df['predicted_value'].isna().any():
        print("Missing predictions found -> skipping (pearson_r = NaN).")
        return None
    
    # Pair forward vs RC on the shared base id (name minus the trailing '_RC'), NOT by
    # row position
    is_rc = predictions_df['seq_name'].str.lower().str.endswith('_rc')
    forward_df = predictions_df[~is_rc].copy()
    rc_df = predictions_df[is_rc].copy()
    rc_df['seq_name'] = rc_df['seq_name'].str.slice(stop=-3)  # strip trailing '_RC'
    
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

    if len(merged) < 1:
        print("No forward/RC pairs -> cannot compute correlation (pearson_r = NaN)")
        return None
    
    correlations = []
    for _, row in merged.iterrows():
        seq_id = row['seq_name']
    
        try:
            fwd = np.asarray(row['predicted_value_fwd'], dtype=float)
            # Reverse the RC track so it lines up position-by-position with the forward
            rc_rev = np.asarray(row['predicted_value_rc'], dtype=float)[::-1]
        except (ValueError, TypeError):
            print(f"Non-numeric track for '{seq_id}' -> skipping (pearson_r = NaN)")
            return None
    
        # pearsonr needs paired 1-D tracks of equal length
        if fwd.ndim != 1 or rc_rev.ndim != 1:
            print(f"Track for '{seq_id}' is not 1-D (forward ndim={fwd.ndim}, "
                  f"RC ndim={rc_rev.ndim}) -> skipping (pearson_r = NaN)")
            return None
        if fwd.shape != rc_rev.shape:
            print(f"Track length mismatch for '{seq_id}' (forward {fwd.shape}, "
                  f"RC {rc_rev.shape}) -> skipping (pearson_r = NaN)")
            return None
        if np.isnan(fwd).any() or np.isnan(rc_rev).any():
            print(f"NA values in track for '{seq_id}' -> skipping (pearson_r = NaN)")
            return None
    
        # Zero variance on either side -> pearsonr is undefined for this pair. The
        # predictions exist, so per the framework convention this pair contributes 0.0
        # rather than poisoning the whole mean with NaN.
        if fwd.std() == 0 or rc_rev.std() == 0:
            correlations.append(0.0)
            continue

        r, _ = pearsonr(fwd, rc_rev)
        correlations.append(0.0 if np.isnan(r) else r)

    mean_r = float(np.mean(correlations))
    print(f"Calculated Mean Pearson r: {mean_r}")
    return mean_r

def calculate_and_save_metrics(saved_predictions_path, output_dir):
    """
    Compute the mean RC-consistency Pearson r across tracks and append one row to the
    summary file evaluation_summary_<EVALUATOR_NAME>.csv (tab-separated, append mode)
    """
    
    try:
        if not os.path.exists(saved_predictions_path):
            print(f"Predictions file not found: {saved_predictions_path}", file=sys.stderr)
            return

        print("Starting evaluation calculation and saving as CSV")
        print(f"Using predictions from: {saved_predictions_path}")
        print(f"Summary will be written to: {output_dir}")

        with open(saved_predictions_path, 'r') as f:
            predictions_file_content = json.load(f)

        # Guard the predictor name so a missing key can't AttributeError on .replace().
        predictor_name_received = predictions_file_content.get("predictor_name") or "UnknownPredictor"
        predictor_name = predictor_name_received.replace(" ", "_").replace("/", "_")

        description = "Consistency track (K562)"
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")

        pearson_r = calculate_pearson_r_track(predictions_file_content)

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