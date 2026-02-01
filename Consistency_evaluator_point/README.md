# Consistency Point Request Evaluator (K562)

A genomic sequence Evaluator that tests consistency of predictions between forward sequences and their reverse complements. This Evaluator requests point chromatin accessibility predictions for K562 cell type (homo_sapiens) and calculates Pearson correlation as a metric of Predictor reliability.

## Overview

This evaluator:
- Sends 100kb Genomic sequences (forward and reverse complement pairs) to Predictors and requests accessibility predictions
- Calculates Pearson correlation coefficient between forward and reverse complement predictions
- Supports both JSON and MessagePack data formats

To get started with running the pre-built container, download from the following link: https://zenodo.org/records/18182295

### Running with Container

```bash
apptainer run --containall \
  --B /path/to/evaluator_data:/evaluator_data \
  --B /path/to/output:/predictions \
  consistency_evaluator_point.sif <predictor_ip> <predictor_port> /predictions
```

**Example:**
```bash
apptainer run --containall \
  --B ./evaluator_data:/evaluator_data \
  --B ./results:/predictions \
  consistency_evaluator_point.sif 192.168.1.100 8000 /predictions
```

If you would like to re-build the container, edit corresponding paths in the .def file and build using the command below. Please make sure to follow the directory stucture shown below. 

### Container Build

```bash
# Build the Apptainer container
apptainer build consistency_evaluator_point.sif evaluator.def
```

## Output Files

The evaluator generates two output files:

### 1. Raw Predictions JSON
`<EVALUATOR_NAME>_predictions_summary_from_<PREDICTOR_NAME>.json`

Contains:
- Predictor name and metadata
- All prediction tasks and their parameters
- Raw predictions for each sequence

### 2. Evaluation Metrics CSV
`<EVALUATOR_NAME>_predictions_summary_from_<PREDICTOR_NAME>.csv`

Contains:
- Evaluator name and description
- Predictor name and timestamp
- Pearson correlation coefficient (r)
- Prediction task metadata

**Example output:**
```tsv
Evaluator_name	Description	Predictor_name	Time_stamp	Metric	Value	Prediction_task(s)_data
Consistency_Point_Request_K562	Consistency point (K562)	MyPredictor	20250125-143022.123456	pearson_r	0.945	[{...}]
```

## API Specification

The evaluator communicates with predictor services via REST API:

### Format Negotiation
```
GET /formats
```
Returns supported request and response formats.

### Prediction Request
```
POST /predict
REQUEST_FORMAT = "application/json"

RESPONSE_FORMAT = "application/json"
```

**Request Body:**
```json
{
  "readout": "point",
  "prediction_tasks": [
    {
      "name": "consistency_K562",
      "type": "accessibility",
      "cell_type": "K562",
      "species": "homo_sapiens"
    }
  ],
  "sequences": {
    "seq_forward": "ATCGATCG...",
    "seq_RC": "CGATCGAT..."
  }
}
```

## Evaluation Metric

**Pearson Correlation Coefficient (r)**

Measures the correlation between forward sequence point predictions and their reverse complement point predictions. 

## Directory Structure

```
Consistency_evaluator_point_Astro/
├── config.py                           # Configuration settings
├── consistency_evaluator_RestAPI.py     # Main evaluator script
├── data_loader.py                       # Input data loading and validation
├── evaluator_content_handler.py         # API communication logic
├── evaluator_metrics_calculator.py      # Metrics computation
├── evaluator.def                        # Apptainer container definition
├── evaluator_data/
│   └── all_consistency_data.csv        # Input sequences
└── README.md                            # This file
```


