#!/usr/bin/env python3
"""
Dadi Optimisation Results Analyser

This script analyses dadi optimisation results to find the best models
based on log-likelihood and AIC for each Pop and Mask combination.

Author: GitHub Copilot
Date: September 2025
"""

import pandas as pd
import numpy as np
import argparse
import sys
from pathlib import Path

def load_data(file_path):
    """Load and clean the dadi optimisation results."""
    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"Loaded {len(df)} rows from {file_path}")
        return df
    except Exception as e:
        print(f"Error loading file: {e}")
        sys.exit(1)

def clean_data(df):
    """Clean the data and handle any formatting issues."""
    # Convert numeric columns
    numeric_cols = ['log-likelihood', 'AIC', 'chi-squared', 'theta']
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Remove rows with missing critical values
    initial_rows = len(df)
    df = df.dropna(subset=['log-likelihood', 'AIC'])
    if len(df) < initial_rows:
        print(f"Removed {initial_rows - len(df)} rows with missing log-likelihood or AIC values")
    
    return df

def extract_pop_mask_combinations(df):
    """Extract unique Pop and Mask combinations."""
    combinations = df[['Pop', 'Mask']].drop_duplicates().sort_values(['Pop', 'Mask'])
    print(f"\nFound {len(combinations)} unique Pop-Mask combinations:")
    for _, row in combinations.iterrows():
        print(f"  - {row['Pop']} | {row['Mask']}")
    return combinations

def find_best_models(df):
    """Find best models by log-likelihood and AIC for each Pop-Mask combination, including second-best."""
    results = []
    
    # Group by Pop and Mask
    grouped = df.groupby(['Pop', 'Mask'])
    
    for (pop, mask), group in grouped:
        # Top 2 runs overall (regardless of model)
        top_runs_by_ll = group.sort_values('log-likelihood', ascending=False)
        top_runs_by_aic = group.sort_values('AIC', ascending=True)
        
        best_run_ll = top_runs_by_ll.iloc[0]
        second_run_ll = top_runs_by_ll.iloc[1] if len(top_runs_by_ll) > 1 else None
        
        best_run_aic = top_runs_by_aic.iloc[0]
        second_run_aic = top_runs_by_aic.iloc[1] if len(top_runs_by_aic) > 1 else None
        
        # Top 2 models (best run of each model type)
        model_best = group.groupby('Model').apply(
            lambda x: x.loc[x['log-likelihood'].idxmax()]
        ).reset_index(drop=True)
        
        model_best_by_ll = model_best.sort_values('log-likelihood', ascending=False)
        model_best_by_aic = model_best.sort_values('AIC', ascending=True)
        
        best_model_ll = model_best_by_ll.iloc[0]
        second_model_ll = model_best_by_ll.iloc[1] if len(model_best_by_ll) > 1 else None
        
        best_model_aic = model_best_by_aic.iloc[0]
        second_model_aic = model_best_by_aic.iloc[1] if len(model_best_by_aic) > 1 else None
        
        # Check if same model wins both criteria
        same_model = (best_model_ll['Model'] == best_model_aic['Model'])
        
        result = {
            'Pop': pop,
            'Mask': mask,
            'N_models_tested': len(model_best),  # Number of unique models
            'N_total_runs': len(group),  # Total number of runs
            'Models_tested': ', '.join(sorted(model_best['Model'].unique())),
            
            # Best RUN by log-likelihood
            'Best_Run_LL_Model': best_run_ll['Model'],
            'Best_Run_LL_Value': best_run_ll['log-likelihood'],
            'Best_Run_LL_AIC': best_run_ll['AIC'],
            'Best_Run_LL_Params': best_run_ll['optimised_params'],
            
            # Second best RUN by log-likelihood
            'Second_Run_LL_Model': second_run_ll['Model'] if second_run_ll is not None else None,
            'Second_Run_LL_Value': second_run_ll['log-likelihood'] if second_run_ll is not None else None,
            'Second_Run_LL_AIC': second_run_ll['AIC'] if second_run_ll is not None else None,
            'Run_LL_Difference': best_run_ll['log-likelihood'] - second_run_ll['log-likelihood'] if second_run_ll is not None else None,
            'Run_Same_Model': best_run_ll['Model'] == second_run_ll['Model'] if second_run_ll is not None else None,
            
            # Best MODEL by log-likelihood
            'Best_Model_LL_Model': best_model_ll['Model'],
            'Best_Model_LL_Value': best_model_ll['log-likelihood'],
            'Best_Model_LL_AIC': best_model_ll['AIC'],
            'Best_Model_LL_Params': best_model_ll['optimised_params'],
            
            # Second best MODEL by log-likelihood
            'Second_Model_LL_Model': second_model_ll['Model'] if second_model_ll is not None else None,
            'Second_Model_LL_Value': second_model_ll['log-likelihood'] if second_model_ll is not None else None,
            'Second_Model_LL_AIC': second_model_ll['AIC'] if second_model_ll is not None else None,
            'Model_LL_Difference': best_model_ll['log-likelihood'] - second_model_ll['log-likelihood'] if second_model_ll is not None else None,
            
            # Best MODEL by AIC
            'Best_Model_AIC_Model': best_model_aic['Model'],
            'Best_Model_AIC_Value': best_model_aic['AIC'],
            'Best_Model_AIC_LL': best_model_aic['log-likelihood'],
            'Best_Model_AIC_Params': best_model_aic['optimised_params'],
            
            # Second best MODEL by AIC
            'Second_Model_AIC_Model': second_model_aic['Model'] if second_model_aic is not None else None,
            'Second_Model_AIC_Value': second_model_aic['AIC'] if second_model_aic is not None else None,
            'Second_Model_AIC_LL': second_model_aic['log-likelihood'] if second_model_aic is not None else None,
            'Model_AIC_Difference': second_model_aic['AIC'] - best_model_aic['AIC'] if second_model_aic is not None else None,
            
            # Summary
            'Same_Best_Model': same_model,
            'Model_LL_AIC_Difference': best_model_ll['AIC'] - best_model_aic['AIC'] if not same_model else 0
        }
        
        results.append(result)
    
    return pd.DataFrame(results)

def calculate_model_summary(df):
    """Calculate summary statistics for each model type."""
    model_stats = []
    
    for model in df['Model'].unique():
        model_data = df[df['Model'] == model]
        
        model_stats.append({
            'Model': model,
            'N_runs': len(model_data),
            'Mean_LL': model_data['log-likelihood'].mean(),
            'Std_LL': model_data['log-likelihood'].std(),
            'Mean_AIC': model_data['AIC'].mean(),
            'Std_AIC': model_data['AIC'].std(),
            'Best_LL': model_data['log-likelihood'].max(),
            'Best_AIC': model_data['AIC'].min()
        })
    
    return pd.DataFrame(model_stats).sort_values('Mean_LL', ascending=False)

def print_summary_report(best_models_df, model_summary_df, title="ALL POPULATIONS"):
    """Print a comprehensive summary report."""
    print("\n" + "="*80)
    print(f"DADI OPTIMISATION RESULTS SUMMARY - {title}")
    print("="*80)
    
    # Overall statistics
    print(f"\nAnalysed {len(best_models_df)} Pop-Mask combinations")
    print(f"Total models tested: {', '.join(sorted(model_summary_df['Model'].unique()))}")
    
    # Best performing models overall
    print(f"\nMODEL PERFORMANCE SUMMARY:")
    print("-" * 50)
    print(model_summary_df.to_string(index=False, float_format='%.2f'))
    
    # Cases where LL and AIC agree/disagree
    same_count = best_models_df['Same_Best_Model'].sum()
    total_count = len(best_models_df)
    print(f"\nMODEL SELECTION AGREEMENT:")
    print("-" * 30)
    print(f"Log-likelihood and AIC agree: {same_count}/{total_count} ({same_count/total_count*100:.1f}%)")
    print(f"Log-likelihood and AIC disagree: {total_count-same_count}/{total_count} ({(total_count-same_count)/total_count*100:.1f}%)")
    
    # Show disagreement cases
    if same_count < total_count:
        print(f"\nCASES WHERE LL AND AIC DISAGREE:")
        print("-" * 40)
        disagreements = best_models_df[~best_models_df['Same_Best_Model']]
        for _, row in disagreements.iterrows():
            print(f"  {row['Pop']} | {row['Mask']}:")
            print(f"    Best by LL: {row['Best_Model_LL_Model']} (LL={row['Best_Model_LL_Value']:.2f}, AIC={row['Best_Model_LL_AIC']:.2f})")
            print(f"    Best by AIC: {row['Best_Model_AIC_Model']} (LL={row['Best_Model_AIC_LL']:.2f}, AIC={row['Best_Model_AIC_Value']:.2f})")
            print(f"    AIC difference: {row['Model_LL_AIC_Difference']:.2f}")
    
    print(f"\nDETAILED RESULTS BY POP-MASK COMBINATION:")
    print("-" * 50)
    
    # Show detailed results for each combination
    for _, row in best_models_df.iterrows():
        print(f"\n{row['Pop']} | {row['Mask']} ({row['N_models_tested']} models, {row['N_total_runs']} total runs)")
        
        # Show top 2 runs
        print(f"  TOP RUNS:")
        print(f"    1st: {row['Best_Run_LL_Model']} (LL: {row['Best_Run_LL_Value']:.2f})")
        if row['Second_Run_LL_Model'] is not None:
            print(f"    2nd: {row['Second_Run_LL_Model']} (LL: {row['Second_Run_LL_Value']:.2f}, diff: {row['Run_LL_Difference']:.2f})")
            consistency = "✓ Same model" if row['Run_Same_Model'] else "⚠ Different models"
            print(f"    Consistency: {consistency}")
        
        # Show top 2 models
        print(f"  TOP MODELS:")
        print(f"    1st: {row['Best_Model_LL_Model']} (LL: {row['Best_Model_LL_Value']:.2f})")
        if row['Second_Model_LL_Model'] is not None:
            print(f"    2nd: {row['Second_Model_LL_Model']} (LL: {row['Second_Model_LL_Value']:.2f}, diff: {row['Model_LL_Difference']:.2f})")
        
        if row['Same_Best_Model']:
            print(f"  ✓ Same model wins by both log-likelihood and AIC")
        else:
            print(f"  ⚠ Different models win by LL ({row['Best_Model_LL_Model']}) vs AIC ({row['Best_Model_AIC_Model']})")

def filter_data_by_criteria(df):
    """Filter data to only include specific models and mask."""
    # Target models and mask
    target_models = ['bottle', 'size_change', 'bottle_neck', 'snm.1d']
    target_mask = 'mid'
    
    # Filter by mask
    filtered_df = df[df['Mask'] == target_mask].copy()
    print(f"Filtered to mask '{target_mask}': {len(filtered_df)} rows")
    
    # Filter by models
    filtered_df = filtered_df[filtered_df['Model'].isin(target_models)].copy()
    print(f"Filtered to target models {target_models}: {len(filtered_df)} rows")
    
    if len(filtered_df) == 0:
        print("Warning: No data remains after filtering!")
        return filtered_df
    
    # Show what we have
    available_models = sorted(filtered_df['Model'].unique())
    print(f"Available models after filtering: {available_models}")
    
    available_pops = sorted(filtered_df['Pop'].unique())
    print(f"Available populations: {len(available_pops)} populations")
    
    return filtered_df

def filter_use_populations(df):
    """Filter data to only include populations with 'USE' in the name."""
    use_pops = df[df['Pop'].str.contains('USE', case=False, na=False)]
    print(f"\nFiltered to {len(use_pops)} rows with 'USE' populations")
    
    if len(use_pops) > 0:
        unique_use_pops = use_pops['Pop'].unique()
        print(f"USE populations found: {len(unique_use_pops)}")
        for pop in sorted(unique_use_pops):
            print(f"  - {pop}")
    
    return use_pops

def save_results(best_models_df, model_summary_df, output_dir, use_best_models_df=None, use_model_summary_df=None):
    """Save results to CSV files."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Save all results
    best_models_file = output_path / 'best_models_by_pop_mask.csv'
    best_models_df.to_csv(best_models_file, index=False)
    print(f"\nAll results - Best models saved to: {best_models_file}")
    
    model_summary_file = output_path / 'model_performance_summary.csv'
    model_summary_df.to_csv(model_summary_file, index=False)
    print(f"All results - Model summary saved to: {model_summary_file}")
    
    # Save USE-only results if available
    if use_best_models_df is not None and len(use_best_models_df) > 0:
        use_best_models_file = output_path / 'best_models_by_pop_mask_USE_only.csv'
        use_best_models_df.to_csv(use_best_models_file, index=False)
        print(f"USE only - Best models saved to: {use_best_models_file}")
        
        if use_model_summary_df is not None and len(use_model_summary_df) > 0:
            use_model_summary_file = output_path / 'model_performance_summary_USE_only.csv'
            use_model_summary_df.to_csv(use_model_summary_file, index=False)
            print(f"USE only - Model summary saved to: {use_model_summary_file}")
    else:
        print("No USE populations found - skipping USE-only results files")

def main():
    parser = argparse.ArgumentParser(description='Analyse dadi optimisation results')
    parser.add_argument('input_file', help='Path to dadi optimisation results file')
    parser.add_argument('-o', '--output', default='.', help='Output directory for CSV files')
    parser.add_argument('--no-save', action='store_true', help='Do not save CSV files')
    
    args = parser.parse_args()
    
    # Load and process data
    print("Loading dadi optimisation results...")
    df = load_data(args.input_file)
    
    print("Cleaning data...")
    df = clean_data(df)
    
    print("Applying filters (mask='mid', models=['bottle', 'size_change', 'bottle_neck', 'snm.1d'])...")
    df_filtered = filter_data_by_criteria(df)
    
    if len(df_filtered) == 0:
        print("No data remaining after filtering. Exiting.")
        sys.exit(1)
    
    print("Extracting Pop-Mask combinations...")
    combinations = extract_pop_mask_combinations(df_filtered)
    
    print("Finding best models for FILTERED populations...")
    best_models_df = find_best_models(df_filtered)
    
    print("Calculating model summary statistics for FILTERED populations...")
    model_summary_df = calculate_model_summary(df_filtered)
    
    # Process USE populations separately
    use_df = filter_use_populations(df_filtered)
    use_best_models_df = None
    use_model_summary_df = None
    
    if len(use_df) > 0:
        print("Finding best models for USE populations only...")
        use_best_models_df = find_best_models(use_df)
        
        print("Calculating model summary statistics for USE populations only...")
        use_model_summary_df = calculate_model_summary(use_df)
    
    # Print summary reports
    print_summary_report(best_models_df, model_summary_df, "FILTERED POPULATIONS (mid mask, target models only)")
    
    if use_best_models_df is not None and len(use_best_models_df) > 0:
        print_summary_report(use_best_models_df, use_model_summary_df, "FILTERED USE POPULATIONS ONLY")
    
    # Save results if requested
    if not args.no_save:
        save_results(best_models_df, model_summary_df, args.output, use_best_models_df, use_model_summary_df)
    
    print(f"\nAnalysis complete!")

if __name__ == "__main__":
    main()