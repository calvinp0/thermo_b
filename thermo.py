import pandas as pd
from typing import List, Callable
import matplotlib.pyplot as plt

def ideal_gas_law(gas_constant: float, list_of_volumes: List, list_of_temperatures: List) -> pd.DataFrame:
    """ A function that calculates the pressure of an ideal gas at a given temperature and volume.

    Args:
    ----
        gas_constant (float): The gas constant - 8.3145 cm^3*MPa/mol*K
        list_of_volumes (List): A list of molar volumes in cm^3mol^-1
        list_of_temperatures (List): A list of temperatures in K

    Returns:
    -------
        pd.DataFrame: A DataFrame with the pressure (MPa) of an ideal gas at a given temperature and molar volume.
    """
    df = pd.DataFrame(columns=list_of_temperatures, index=list_of_volumes)
    df = df.apply(lambda x: (float(x.name) * gas_constant) / x.index)
    return df


# Plot Function for Task 2 - Ideal Gas Law

def plot_ideal_gas_law(ideal_gas_law_df: pd.DataFrame) -> None:
    """ A function that plots the pressure of an ideal gas at a given temperature and volume provided by a DataFrame.

    Args:
    ----
        ideal_gas_law_df (pd.DataFrame): A DataFrame with the pressure (MPa) of an ideal gas at given temperatures and molar volumes.
    """

    for column in ideal_gas_law_df:
        plt.plot(ideal_gas_law_df.index, ideal_gas_law_df[column], label=column)
    plt.ylim(0, 15)
        

# Peng-Robinson EOS Function for Task 3

def peng_robinson_eos(gas_constant: float, list_of_volumes: List, list_of_temperatures: List, bc: float, ac: float, alpha: Callable) -> pd.DataFrame:
    """ A function that calculates the pressure of a gas at a given temperature and molar volume using the Peng-Robinson EOS.

    Args:
        gas_constant (float): 8.3145 cm^3*MPa/mol*K
        list_of_volumes (List): List of molar volumes in cm^3mol^-1
        list_of_temperatures (List): List of temperatures in K
        bc (float): Value of b at the critical temperature
        ac (float): Value of c at the critical temperature
        alpha (Callable): Alpha function is a function of temperature that is used to calculate the acentric factor of the gas

    Returns:
        pd.DataFrame: A DataFrame with the pressure (MPa) of a gas at a given temperature and molar volume using the Peng-Robinson EOS.
    """
    df = pd.DataFrame(columns=list_of_temperatures, index=list_of_volumes)
    temperature_aT = {T: (ac * alpha(T)) for T in list_of_temperatures}
    df = df.apply(lambda x: ((gas_constant * float(x.name)) / (x.index - bc)) - (temperature_aT[float(x.name)] / ((x.index * (x.index + bc) + bc * (x.index - bc)))), axis=0)
    return df

# Plot Function for Task 4 - Peng-Robinson EOS

def plot_peng_robinson_eos(peng_robinson_eos_df: pd.DataFrame) -> None:
    """ A function that plots the pressure of a gas at a given temperature and molar volume provided by a DataFrame.

    Args:
        peng_robinson_eos_df (pd.DataFrame): A DataFrame with the pressure (MPa) of a gas at given temperatures and molar volumes using the Peng-Robinson EOS.
    """
    
    for column in peng_robinson_eos_df:
        plt.plot(peng_robinson_eos_df.index, peng_robinson_eos_df[column], label=column)

# PLot Function for Task 6

def plot_example_ideal_gas_law_and_peng_robinson_eos() -> None:
    """A function that plots the pressure of an ideal gas and a gas using the Peng-Robinson EOS at a given temperature and molar volume provided by an example DataFrame.
    """
    example = pd.read_csv('src/preos_df_example.csv', index_col=0)
    plot_peng_robinson_eos(example)
    plt.axhline(y=4.6, color='red', linestyle='--', label='$P^{0}$ at 284.0 K (Lit)')
    plt.xlabel('Volume (m^3)')
    plt.ylabel('Pressure (MPa)')
    plt.ylim(0, 15)
    plt.title('Peng-Robinson EOS - CO2 Example')
    plt.legend()
    plt.show()


def get_chemical_values(chemical_df: pd.DataFrame, molecules: List, column: str) -> List:
    """ A function that returns the values of a column for a list of molecules.

    Args:
        chemical_df (pd.DataFrame): A DataFrame with the chemical values of molecules
        molecules (List): A list of molecules
        column (str): The column of the DataFrame to be returned

    Returns:
        List: A list of values of the column for the molecules
    """
    if column not in chemical_df.columns:
        raise ValueError('Column not in DataFrame. Did you input the correct column name?')
    elif molecules is None:
        return ValueError('No molecules were inputted. Please input a list of molecules.')
    elif not isinstance(molecules, list):
        raise TypeError('Molecules must be a list.')
    elif [molecule in chemical_df['Element'].unique() for molecule in molecules] != [True, True]:
        raise ValueError('Molecules not in DataFrame. Did you input the correct molecules?')
    elif chemical_df is None:
        raise ValueError('No DataFrame was inputted. Please input the required DataFrame.')
    else:
        return chemical_df.loc[chemical_df['Element'].isin(molecules), column].values.tolist()
    

# Calculate ai for Task 7.1.1
def calculate_ai(chemical_df: pd.DataFrame, molecules: List, temperature: float) -> List:
    """ A function that calculates the ai values for a list of molecules at a given temperature.

    Args:
        chemical_df (pd.DataFrame): A DataFrame with the chemical values of molecules
        molecules (List): Molecules to calculate ai values for
        temperature (float): Temperature in K

    Returns:
        List: A list of ai values for the molecules
    """
    if temperature <= 0:
        raise ValueError('Temperature must be greater than 0.')
    elif not isinstance(temperature, float):
        raise TypeError('Temperature must be a float.')
    else:
        #  0.45724 * R**2 * df.loc['Tc_K']**2 / df.loc['Pc_MPa'] * (1 + df.loc['ki'] * (1 - (T/df.loc['Tc_K'])**0.5))**2
        return 0.45724 * (8.3145**2) * chemical_df.loc['Tc_K']**2/chemical_df.loc['Pc_MPa'] * (1 + chemical_df.loc['ki'] * (1 - (temperature/chemical_df.loc['Tc_K'] )**0.5))**2


# Task 7.2

def calculate_df_aij(molecules: List, df: pd.DataFrame, df_interaction: pd.DataFrame, df_aij: pd.DataFrame):
    for i in molecules:
        for j in molecules:
            df_aij.loc[i, j] = (df.loc['ai', i] * df.loc['ai', j])**0.5 * (1 - df_interaction.loc[i, j])
    return df_aij