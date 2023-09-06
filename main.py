import streamlit as st
import pandas as pd
import re
import plotly.express as px
import numpy as np
from scipy.optimize import curve_fit

uploaded_file = st.file_uploader("upload feature matrix", type="tsv")
if uploaded_file:
    dataframe = pd.read_csv(uploaded_file)
else:
    st.info("💡 No input file specified. Using example data instead.")

df = pd.read_csv("metabolite_intensity_data.tsv", sep="\t", index_col="metabolite")

metabolite = st.selectbox("metabolite", df.index)

def is_standard(s):
    if "mM" in s or "µM" in s or "nM" in s:
        return True
    else:
        return False

standards = st.multiselect("standards", df.columns, [c for c in df.columns if is_standard(c)])

samples = st.multiselect("samples", df.columns, [c for c in df.columns if c not in standards])

conc_unit = st.selectbox("concentration unit", ["mM", "µM", "nM"])

def get_conc_value(s, u):
    factors = {"mM": 1000, "µM": 1000000, "nM": 1000000000}
    matches_num = re.findall(r"\d+[.,]*\d+|0", s)
    matches_unit = re.findall(r"[mµn]M", s)
    if not matches_num or not matches_unit:
        return s
    n = float(matches_num[0].replace(",", "."))
    n /= factors[matches_unit[0]]
    n *= factors[u]
    if n.is_integer():
        n = int(n)
    return str(n)

s = pd.DataFrame({"Sample Name": standards, conc_unit: [get_conc_value(s, conc_unit) for s in standards]}).set_index("Sample Name")

edited = st.data_editor(s, use_container_width=True, disabled=["Sample Name"])

faulty = []
for x in edited[conc_unit]:
    try:
        int_value = float(x)
    except ValueError:
        faulty.append(x)
if faulty:
    st.info("Invalid concentration values found: ")
    for x in faulty:
        st.warning(x)

else:
    df = df.rename(columns={old: new for old, new in zip(edited.index, edited[conc_unit])})
    
    data = df.loc[metabolite, edited[conc_unit]]
    data.index = pd.to_numeric(data.index, errors='coerce')
    data = data.sort_index()

    concs = edited[conc_unit].astype(int).sort_values()
    linear_region = st.select_slider("select a region where standard intensities are linear", concs, (concs.min(), concs.max()))

    data_linear = data.loc[[i for i in data.index if i >= linear_region[0] and i <= linear_region[1]]]

    data_linear = pd.DataFrame(data_linear)

    # Define the linear function for regression
    def linear(x, m, c):
        return m * x + c
    
    # Predict concentration based on intensity
    def predict_linear(x, m, c):
        return (x - c) / m

    # Fit the linear regression line to the data
    popt_linear, _ = curve_fit(linear, data_linear.index, data_linear[metabolite])

    # Generate x-values for the regression lines
    x_values = np.linspace(data_linear.index.min(), data_linear.index.max(), 100)

    # Calculate y-values for the linear regression line
    y_values_linear = linear(x_values, *popt_linear)

    names = [str(concs.loc[n]) if n in concs.index else n for n in samples]
    predicted_values = pd.DataFrame({"Sample": samples, "Intensity": [df.loc[metabolite, s] for s in names]})
    predicted_values[f"Predicted Concentration ({conc_unit})"] = predict_linear(predicted_values["Intensity"], *popt_linear)
    predicted_values = predicted_values.set_index("Sample")
    predicted_values[f"Predicted Concentration ({conc_unit})"] = predicted_values[f"Predicted Concentration ({conc_unit})"].apply(lambda x: int(x) if x > 10 else x)
    predicted_values[f"Predicted Concentration ({conc_unit})"] = predicted_values[f"Predicted Concentration ({conc_unit})"].apply(lambda x: 0 if x < 0 else x)


    # Create a Plotly line plot
    fig = px.line(data_linear, x=data_linear.index, y=metabolite, markers=True, line_shape='linear')
    fig.update_traces(name="standard intensities", showlegend=True)

    # Add the linear regression line to the plot
    fig.add_scatter(x=x_values, y=y_values_linear, mode='lines', name='linear regression', line=dict(color='green'))

    # Plot data points which have been excluded (not linear)
    data_outliers = data.loc[[i for i in data.index if i not in data_linear.index]]
    data_outliers = pd.DataFrame(data_outliers)

    fig.add_scatter(x=data_outliers.index, y=data_outliers[metabolite], mode='markers',
                    marker=dict(symbol='x', size=8, color='red'), line=dict(color='red'), name="excluded values")

    fig.update_layout( xaxis_title=f"concentration {conc_unit}", height=600)
    # Show the plot
    st.plotly_chart(fig, use_container_width=True)

    st.dataframe(predicted_values, use_container_width=True)