# clpm-sleep-health-longitudinal
# Exploring Bidirectional Relationship of Longitudinal Data 
# Using Cross-Lagged Panel Models

## Overview
This project demonstrates the bidirectional longitudinal 
relationship between sleep quality and general health (GHQ-12) 
using five waves of the UK Household Longitudinal Study (UKHLS), 
spanning 2009 to 2021 (N = 12,875).

## Methods
- Cross-Lagged Panel Model (CLPM) using R (lavaan)
- Maximum Likelihood estimation with FIML for missing data
- Constrained vs unconstrained model comparison via LR test
- Multigroup moderation by gender, employment, education
- Descriptive visualisations in Python (matplotlib, seaborn)

## Key Findings
- Bidirectional relationship confirmed at all four wave transitions
- GHQ → Sleep effect strengthened over time (β = 0.055 to 0.131)
- Moderation confirmed for gender, employment, and education

## Software
- R 4.5.2 — lavaan, haven, tidyverse, ggplot2
- Python 3.12 — pandas, numpy, matplotlib, seaborn, statsmodels

## Data Source
Understanding Society: UK Household Longitudinal Study (UKHLS)
Study Number 6614, UK Data Service
Note: Raw data files are not included in this repository 
in accordance with UK Data Service End User Licence terms.
