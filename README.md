## Introduction

Our final project mainly contains two parts: a traditional SEIR model with theoretical analysis which can give users an understanding of how parameters and initial value will affect the solution of an ordinal derivative equation through interaction with Shiny App. A complicated Age Stratification SEIR model that considers more practical situations simulates the spread of COVID-19.

For the first part, we give a detailed theoretical analysis of the basic SEIR model: Model definition, Meaning of parameters, feasible region, and Equilibrium point. Then we use the basic reproduction number $\mathcal{R}_{0}$ to show the situation of local stability and global stability. In the shiny app interface, users can first upload a dataset about disease spread, then we use the least square method to estimate the parameters in the model. After that, we will simulate the spread of disease through the SEIR model and draw several plots of trajectories and phase planes. By combining theoretical analysis and numerical simulation, users can have a deeper understanding of the SEIR model which is an ordinary differential equation. For example, we provide toy data, the estimated parameters make $\mathcal{R}_{0}$ larger than 1. Using the theoretical analysis, we know that the unique endemic equilibrium of the system is globally asymptotically stable in the interior of the feasible region. However, the plot in the simulation part doesn't show a stable trend with the time from 0 to 1000. Then the theory can guide numerical simulation, as users drag the time slider, it finally converges in the trajectories and phase planes.

For the second part, we use a proposed cutting-aged SEIR model variant [1]: SEIcIscR Model. It is an age-structured SEIR model separating the I Compartment into to Ic and Isc Compartment so that it can define different contact rates between different age-group and consider contributions of asymptomatic and sub-clinical cases. Users can also change the parameters and see the trajectories of each compartment with a specific age group.

## App Structure

### Packages

`deSolve` is required to solve differential equations.\
`tidyverse`, `ggplot2` is used to perform data manipulation and visualization\
`tableone`, `kableExtra` is used to output table.\
`ggpubr`, `ggsci`, `Polychrome` is used to polish plots.\
`plotly` is used to construct interactable plots.\
`shiny`, `shinythemes`, `shinydashboard`, `dashboardthemes`, `shinyjs` is used to construct shiny app.

### Models

To simulate the transmission among S, E, I, and R compartments, two models are built-in functions. One only includes the basic SEIR model to illustrate its characteristics, while the other added age stratification to further present the distinct trends of different compartments between groups.

1.  Basic SEIR model

    This model mainly describes the population movement among different compartments. It receives `time` (maximum simulation period), `Y` (initial state), and `params` as inputs, returning a matrix that records the population of each compartment across time. Parameters are consist of $\epsilon$ (latency rate), $\lambda$ (contact rate), $\gamma$ (recovery rate) and $\mu$ (ground death rate). Users can freely set their values to explore the conditions that lead to the equilibrium of the model by checking the phase plane or convergence by looking into the trajectories. The default value of $\lambda$ is also provided by estimation via the self-defined `MSE` function and `optimize` function when users load empirical data.

    ``` r
    SEIR <- function(time, Y, params)
    ```

2.  Age-structured SEIR model

    At this stage, age stratification is taken into account. That is to say, contact rates among groups may vary. Therefore, contact rate $\lambda$ is replaced by the product of transmission rate $\beta$ and contact matrices $C$. $C_{i,j}$ denotes contacts of age group $i$ made by age group $j$. The age-structured SEIR function receives the same input as the basic model does, while its `params` may need some modification: $\alpha$ (proportion of infectiousness), $\beta$ (transmission rate), $d_L$ (average incubation period) and $d_I$ (average duration of infection) can be adjusted by the users. There are also other parameters that are collected from other papers and cannot be set: $C$ (contact matrices), $\rho$ (a proportion that infected cases are clinical), $\kappa$ (infectious rate), and $\gamma$ (recovery rate).

    The function returns the population of each compartment across time and can be applied to further compare the dynamic process between different age groups.

    ``` r
    SEIR_age <- function(time, Y, params)
    ```
    
### Esitimation

1.  MSE

    For the optimization of the baseline model, MSE function is defined as the cost function, and `optimize` is called to generate the $\hat{\beta}$. which minimizes the MSE.

    ``` r
    mse <- function(lambda, data, func, times, 
                    initial_state, params)
    ```

### Shiny App

The app is constructed with `dashBoard` in 2 parts: main body and side bar. `tabPanel` is used in main body to construct 3 tabs for different outcome outputs. `ConditionPanle` is used in sidebar so that the UI elements can change with different main body tabs. 

Parameters for the first tab's output is stored in `params_reactive` and parameters for the second tab output is stored in `params_reactive_2.1` and `params_reactive_2.2` corresponding to scenario 1 and 2. The corresponding simulation outcomes are stored in `result_reactive`, `result_reactive_2.1` and `result_reactive` . `updateSliderInput` is used to update the parameters setting in scenario 1 as user change the parameters in the first tab. All there parameters in reactive to the selected model so that they can all be used regardless of selected model.

Eventhough only one slider are provided for time input, the simulation runs in 2 different time scale, one for the real time input value, the other multiply the value by 100. A larger time scale will seriously reduce the speed of the simulation and plots generation, but a lower time scale will be insufficient to see the equilibrium. Results for the interaction plots will only used the real time value as large time value will lag the plot, and result for non-interaction plots and table used the augmented time value for user to see the equilibrium. All simulation of SEIcIscR model is runned under true time value as the simulation for this model is complicated and time consuming.

When plotting the interactable plots, the data need to be accumulated so that user can see the evolution of the lines. Instead of accumulating the real time value, we use time chunks (calculated by divide the time by 10 and remove the decimals) to accumulate the data as it will the frame of the plots and reduce time to construct the interactable plots.

## Instruction for Shiny App

Our shiny app will automatically run the simulation based on the default parameters. Users can choose different preset parameters in **Choose Model parameters**. Users can also choose to manually adjust the different parameters and the app will run the simulation based on updated parameters.

The app also allows users to upload data and it will estimate the **Contact Rate** of the data and run the simulation based on the estimated parameters.

Users can choose different models to do the simulation. The app's default simulation model is the SEIR model, and the user can change to the SEIcIscR model with age saturated structure in **Choose Model**, and the parameters slider will change based on the type of input data of different models.

In the first tab, the user can change the time for the simulation and see the simulation outcome. When the SEIR model is selected, the user will be able to see the changing of different compartments of the model in the first plot. Because the time sometimes is too small for the model equilibrium, the second plot offers a larger time scale (times the time with 100) for the user to see the equilibrium. If the SEIcIsc model is selected, the plot will be presented in 5 chambers, representing 5 different compartments of the model. Both the first Plot in the SEIR model and the plot in the SEIcIscR model are intractable plots, in which the user can choose to hide or show lines in the plot and see the model outcome at different times.

In the third tab, we provide a detailed document about the mathematical expression of the model and theoretical analysis depending on which model the user chooses. Users can gain a basic understanding of differential equations and the derivation of some properties which is instruction and interpretation when they use the interaction function of our Shiny App. For example, for the SEIR model, we provide the model definition, meaning of parameters, feasible region, and equilibrium point. And some properties related to basic reproduction number $\mathcal{R}_{0}$, especially local stability and global stability. Through this knowledge, users can learn how to adjust the coefficients and understand outcomes. As shown above, the user understands that the SEIR model has global stability, and when $\mathcal{R}_{0} \ge 1$, it will be stable at endemic equilibrium $P^{*}=(S^{*},E^{*},I^{*})$. This can be used as a guide to prompt them to slide the time bar on the second plot.

## Reference

Bjørnstad O N, Shea K, Krzywinski M, et al. The SEIRS model for infectious disease dynamics[J]. Nature methods, 2020, 17(6): 557-559.

Bjørnstad O N, Shea K, Krzywinski M, et al. Modeling infectious epidemics[J]. Nature methods, 2020, 17(5): 455-456.

Shea K, Bjørnstad O N, Krzywinski M, et al. Uncertainty and the management of epidemics[J]. Nature methods, 2020, 17(9): 867-869.
