import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import os
import scipy.io
from tools.GSSA_countries import GSSA_country_df
from tools.GSSA_1_agent import GSSA_ShowcaseResult


def show_values_on_bars(axs, h_v="v", space=0.4):
    '''
    A small function that shows the exact integral values of the barplots.
    ------
    Arguments:
        axs(matplotlib.axes): Target subplot.
        h_v(str): Whether the barplot is horizontal or vertical. Default is "v".
        space(float): The space between the tex and the top edge of the bar.

    '''
    def _show_on_single_plot(ax):
        if h_v == "v":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() / 2
                _y = p.get_y() + p.get_height()
                value = int(p.get_height())
                ax.text(_x, _y, value, ha="center") 
        elif h_v == "h":
            for p in ax.patches:
                _x = p.get_x() + p.get_width() + float(space)
                _y = p.get_y() + p.get_height()
                value = int(p.get_width())
                ax.text(_x, _y, value, ha="left")

    if isinstance(axs, np.ndarray):
        for idx, ax in np.ndenumerate(axs):
            _show_on_single_plot(ax)
    else:
        _show_on_single_plot(axs)


def Figure1():
    '''
    This function showcase the ergodic sets
    '''
    # importing data
    real_path = os.path.join(os.getcwd(), "data\\")
    prep1 = scipy.io.loadmat(real_path+r"country_k.mat").get("k")
    prep2 = scipy.io.loadmat(real_path+r"aT20200N10.mat").get("a20200")
    df1 = pd.DataFrame(prep1)
    df2 = pd.DataFrame(prep2)
    # cleaning data
    df1 = df1.iloc[:2000, :1]
    df2 = df2.iloc[:, :1]
    df2 = df2.head(2000)
    result = pd.concat([df1[0], df2[0]], axis=1)
    result.columns = ["Capital", "Productivity"]
    # ploting
    Fig1 = sns.jointplot("Capital", "Productivity", data=result, kind='hex')
    Fig1.fig.suptitle("Ergodic set with closed-formed solution", fontsize=16)
    Fig1.fig.subplots_adjust(top=0.95)
    return


def Figure3():
    '''
    This function showcase the normal polynominal degree vs hermite ones
    '''
    # Define the x for plotting
    x = np.linspace(-2, 2, 100000)

    # The ordinary polynominal degree
    y = x

    # setting the axes at the centre
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plt.subplots_adjust(wspace=0.5)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.title.set_text("Ordinary polynomials $O_{m}(x)$")

    # Plot the functions
    for i in range(6):
        plt.plot(x, y**i, label='$O_{%i}(x)$' %i)

    # Add legend
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)

    # Hermite polynominal degree
    y0 = x**0
    y1 = x**1
    y2 = x**2-1
    y3 = x**3-3*x
    y4 = x**4-6*x**2+3
    y5 = x**5-10*x**3+15*x

    # Setting the axes again at the centre
    ax = fig.add_subplot(1, 2, 2)
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.title.set_text("Hermite polynomials $H_{m}(x)$")

    # Plot the functions
    plt.plot(x, y0, label='$H_{0}(x)$')
    plt.plot(x, y1, label='$H_{1}(x)$')
    plt.plot(x, y2, label='$H_{2}(x)$')
    plt.plot(x, y3, label='$H_{3}(x)$')
    plt.plot(x, y4, label='$H_{4}(x)$')
    plt.plot(x, y5, label='$H_{5}(x)$')

    # Add legend
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)

    # Show the plot
    plt.show()
    return


def Figure2(showcase_result):
    '''
    Benchmarking results for showcase GSSA
    ------
    Arguments:
        showcase_result(pandas DF): The result from GSSA_showcase().
    '''
    colors = ["lightcoral", "indianred", "brown", "firebrick", "maroon"]
    fig = plt.figure(figsize=(10, 5))
    plt.subplots_adjust(wspace=0.5)

    # Subplot 1: Total time
    ax = fig.add_subplot(1, 3, 1)
    plt.bar("Polynomial Degree", "Total Time", width=0.4, data=showcase_result, color=colors)
    plt.ylabel("Total Time(sec)")
    plt.xlabel("Polynomial Degree")

    # Subplot 2: Mean Error
    ax = fig.add_subplot(1, 3, 2)
    plt.bar("Polynomial Degree", "Mean Error", data=showcase_result, color=colors)
    plt.gca().invert_yaxis()
    plt.xlabel("Polynomial Degree")
    plt.ylabel("Mean Error($\log_{10}$)")
    ax = fig.add_subplot(1, 3, 3)
    plt.bar("Polynomial Degree", "Maximum Error", data=showcase_result, color=colors)
    plt.gca().invert_yaxis()
    plt.xlabel("Polynomial Degree")
    plt.ylabel("Maximum Error($\log_{10}$)")

    # build legend
    labels = ["Polynomial Degree "+str(i) for i in range(1, 6)]
    handles = [plt.Rectangle((0, 0), 1,1, color=colors[i]) for i in range(5)]
    plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    plt.suptitle("GSSA Showcase Benchmark", fontsize=16)
    plt.show()
    return


def Figure4():
    '''
    Showcase of inaccuracy of Monte-Carlo intergration method.
    '''
    # Change the IM in showcase
    IM_MonteCarlo = GSSA_ShowcaseResult(IM=0)
    IM_Q1 = GSSA_ShowcaseResult(IM=1)
    IM_Q10 = GSSA_ShowcaseResult(IM=10)
    # Data cleaning
    IM_MonteCarlo['Integration Method'] = "Monte-Carlo Integration"
    IM_Q1['Integration Method'] = "1-Node Gauss-Hermite Integration"
    IM_Q10['Integration Method'] = "10-Node Gauss-Hermite Integration"
    IM_final = pd.concat([IM_MonteCarlo, IM_Q1, IM_Q10])
    # Plot
    fig = plt.figure(figsize=(14, 5))
    fig.suptitle("Figure: Comparision between different intergration methods", fontsize=16)
    plt.subplots_adjust(wspace=0.5)
    ax1 = fig.add_subplot(1, 3, 1)
    g = sns.barplot(x="Polynomial Degree", y="Total Time", hue="Integration Method", data=IM_final, ax=ax1)
    show_values_on_bars(ax1)
    ax1.legend().remove()
    ax2 = fig.add_subplot(1, 3, 2)
    g = sns.barplot(x="Polynomial Degree", y="Mean Error", hue="Integration Method", data=IM_final, ax=ax2)
    g.invert_yaxis()
    g.set(ylabel="Mean Error($\log_{10}$)")
    ax2.legend().remove()
    ax3 = fig.add_subplot(1, 3, 3)
    g = sns.barplot(x="Polynomial Degree", y="Maximum Error", hue="Integration Method", data=IM_final, ax=ax3)
    g.set(ylabel="Maximum Error($\log_{10}$)")
    g.invert_yaxis()
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    plt.show()
    return


def LS_Figure(x, y):
    '''
    Plots the comparison between LS methods.
    -----
    Arguments:
        x(pandas Dataframe): Result1 from the function Result_agent().
        y(pandas Dataframe): Result2 from the function Result_agent().
    '''
    # Get the results
    result1 = x
    result2 = y
    # Some cleaning
    result1.iloc[2:5,3]=0
    result1.iloc[8:10,3]=0
    # Preparation for plotting
    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    plt.subplots_adjust(wspace=0.5)
    fig.suptitle("Figure: LS Methods Comparison", fontsize=16)

    # Plotting
    g = sns.barplot(x="Polynomial Degree", y="Total Time", hue="Method", data=result1, ax=axes[0, 0])
    show_values_on_bars(axes[0, 0])
    axes[0, 0].legend().remove()
    g.set(ylabel="Total Time(sec)")
    g = sns.barplot(x="Polynomial Degree", y="Mean Error", hue="Method", data=result1, ax=axes[0, 1])
    g.invert_yaxis()
    axes[0, 1].legend().remove()
    g.set(ylabel="Mean Error($\log_{10}$)")
    g = sns.barplot(x="Polynomial Degree", y="Max Error", hue="Method", data=result1, ax=axes[0, 2])
    g.invert_yaxis()
    axes[0, 2].legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    g.set(ylabel="Maximum Error($\log_{10}$)")
    g = sns.barplot(x="Polynomial Degree", y="Total Time", hue="Method", data=result2, ax=axes[1, 0])
    show_values_on_bars(axes[1, 0])
    axes[1, 0].legend().remove()
    g.set(ylabel="Total Time(sec)")
    g = sns.barplot(x="Polynomial Degree", y="Mean Error", hue="Method", data=result2, ax=axes[1, 1])
    g.invert_yaxis()
    axes[1, 1].legend().remove()
    g.set(ylabel="Mean Error($\log_{10}$)")
    g = sns.barplot(x="Polynomial Degree", y="Max Error", hue="Method", data=result2, ax=axes[1, 2])
    g.invert_yaxis()
    g.set(ylabel="Maximum Error($\log_{10}$)")
    axes[1, 2].legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    return


def LAD_figure(x):
    '''
    Plots the comparison between LAD methods.
    -----
    Arguments:
        x(pandas Dataframe): Result3 from the function Result_agent().
    '''
    result3 = x
    # Data minor cleaning
    result3.iloc[4, 3] = 0
    result3.iloc[9, 3] = 0

    # Plot
    fig = plt.figure(figsize=(14, 5))
    fig.suptitle("Figure: LAD-PP Methods Comparison", fontsize=16)
    plt.subplots_adjust(wspace=0.5)
    ax1 = fig.add_subplot(1, 3, 1)
    g = sns.barplot(x="Polynomial Degree", y="Total Time", hue="Method", data=result3, ax=ax1)
    g.set(ylabel="Total Time(sec)")
    show_values_on_bars(ax1)
    ax1.legend().remove()
    ax2 = fig.add_subplot(1, 3, 2)
    g = sns.barplot(x="Polynomial Degree", y="Mean Error", hue="Method", data=result3, ax=ax2)
    g.invert_yaxis()
    g.set(ylabel="Mean Error($\log_{10}$)")
    ax2.legend().remove()
    ax3 = fig.add_subplot(1, 3, 3)
    g = sns.barplot(x="Polynomial Degree", y="Max Error", hue="Method", data=result3, ax=ax3)
    g.invert_yaxis()
    g.set(ylabel="Maximum Error($\log_{10}$)")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    return


def country_Figure1():
    '''
    Produce figure for 1 and 2 countries with polynomial degree 1 to 5.
    '''
    countries_1 = GSSA_country_df(N=1, Cache=True)
    countries_2 = GSSA_country_df(N=2, Cache=True)
    countries_12 = pd.concat([countries_1, countries_2])
    # Make sure data in the correct type
    countries_12["Number of countries"] = countries_12["Number of countries"].astype(str)
    # Plotting
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Figure: Comparison between 1 and 2 countries", fontsize=16)
    plt.subplots_adjust(wspace=0.4)
    g = sns.barplot(x="Polynomial Degree", y="Total Time", hue="Number of countries", data=countries_12, palette="rocket_r", ax = ax0)
    show_values_on_bars(ax0)
    ax0.legend().remove()
    g = sns.barplot(x="Polynomial Degree", y="Mean Error", hue="Number of countries", data=countries_12, palette="rocket_r", ax = ax1)
    g.invert_yaxis()
    g.set(ylabel="Mean Error($\log_{10}$)")
    ax1.legend().remove()
    g = sns.barplot(x="Polynomial Degree", y="Max Error", hue="Number of countries", data=countries_12, palette="rocket_r", ax = ax2)
    g.invert_yaxis()
    g.set(ylabel="Maximum Error($\log_{10}$)")
    ax2.legend(loc='center left', title="Number of Countries", bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    plt.show()
    return


def country_Figure2():
    '''
    Produce figure for 10 to 100 countries.
    '''
    df = GSSA_country_df(N=10, Cache=True)
    for i in range(20, 110, 10):
        df_temp = GSSA_country_df(N=i, Cache=True)
        df = df.append(df_temp)
    # plot
    fig = plt.figure(figsize=(15, 5))
    fig.suptitle("Figure: Comparison between 10 to 100 countries", fontsize=16)
    plt.subplots_adjust(wspace=0.5)
    ax1 = fig.add_subplot(1, 3, 1)
    g = sns.barplot(x="Number of countries", y="Total Time", data=df, ax=ax1)
    show_values_on_bars(ax1)
    g.set(xlabel="Number of Countries", ylabel="Total Time (sec)")
    ax2 = fig.add_subplot(1, 3, 2)
    g = sns.barplot(x="Number of countries", y="Mean Error", data=df, ax=ax2)
    g.invert_yaxis()
    g.set(xlabel="Number of Countries", ylabel="Mean Error($\log_{10}$)")
    ax3 = fig.add_subplot(1, 3, 3)
    g = sns.barplot(x="Number of countries", y="Max Error", data=df, ax=ax3)
    g.invert_yaxis()
    g.set(xlabel="Number of Countries", ylabel="Maximum Error($\log_{10}$)")
    plt.show()
    return
