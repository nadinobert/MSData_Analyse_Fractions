# This function has a 2 dim df as input
# Output are the numbers of column 1 and column 2 when ever the number in column 2 changes
def get_flowrate_changes(df):
    # Initialize empty lists to store time points and corresponding flow rates
    result = []

    # Iterate through the rows of the DataFrame
    for i in range(0, len(df)):
        # Check if the flow rate at the current row is different from the previous row
        if i == 0 or df.iloc[i]['flow_rate'] != df.iloc[i - 1]['flow_rate']:
            result.append([df.iloc[i]['time_point'], df.iloc[i]['flow_rate']])

    return result


# This function converts timepoints [min] to elution volume according to the flowrate
def time_to_elution_volume(time, flowrate_list):
    elution_volume = 0.0
    for index, [time_changed, flow] in enumerate(flowrate_list):
        if index == len(flowrate_list) - 1:
            elution_volume += (time - time_changed) * flow
            break
        if time <= time_changed:
            break
        next_time_changed = flowrate_list[index + 1][0]
        if time >= next_time_changed:
            elution_volume += (next_time_changed - time_changed) * flow
        else:
            elution_volume += (time - time_changed) * flow
            break
    return elution_volume

