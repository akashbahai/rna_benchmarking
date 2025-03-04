import os
import pandas as pd

# Directory containing CSV files
csv_directory = './'

# Output CSV file
output_csv = 'interaction_all.csv'

# List all CSV files in the directory that start with 'inf'
csv_files = [f for f in os.listdir(csv_directory) if f.endswith('.csv') and f.startswith('inf')]

# Create an empty DataFrame to store the results
result_df = pd.DataFrame(columns=['Target', 'df', 'drfold', 'rf', 'rhofold', 'rnac', 'trr', '3drna'])

# Process each CSV file
for csv_file in csv_files:
    # Read the CSV file into a DataFrame, explicitly setting 'target' as the index
    df = pd.read_csv(os.path.join(csv_directory, csv_file), index_col='target')

    # Extract target name from the filename
    target_name = csv_file.split('_')[1].split('.')[0]

    # Create a row for the result DataFrame
    row = pd.DataFrame({
        'Target': [target_name],
        'df': ['x'],
        'drfold': ['x'],
        'rf': ['x'],
        'rhofold': ['x'],
        'rnac': ['x'],
        'trr': ['x'],
        '3drna': ['x']
    })

    # Extract inf_all values and update the row if conditions are met
    if df['fn'].str.contains('_df').any():
        row['df'] = df.loc[df['fn'].str.contains('_df'), 'inf_all'].values[0]

    if df['fn'].str.contains('_drfold').any():
        row['drfold'] = df.loc[df['fn'].str.contains('_drfold'), 'inf_all'].values[0]

    if df['fn'].str.contains('_rf').any():
        row['rf'] = df.loc[df['fn'].str.contains('_rf'), 'inf_all'].values[0]

    if df['fn'].str.contains('_rhofold').any():
        row['rhofold'] = df.loc[df['fn'].str.contains('_rhofold'), 'inf_all'].values[0]

    if df['fn'].str.contains('_rnac').any():
        row['rnac'] = df.loc[df['fn'].str.contains('_rnac'), 'inf_all'].values[0]

    if df['fn'].str.contains('_trr').any():
        row['trr'] = df.loc[df['fn'].str.contains('_trr'), 'inf_all'].values[0]

    if df['fn'].str.contains('_3drna').any():
        row['3drna'] = df.loc[df['fn'].str.contains('_3drna'), 'inf_all'].values[0]

    # Append the row to the result DataFrame
    result_df = pd.concat([result_df, row], ignore_index=True)

# Save the result DataFrame to a new CSV file

result_df = result_df.sort_values(by='Target')
result_df.to_csv(output_csv, index=False)
