import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt


def Figure1():
    '''
    This function showcase the ergodic sets
    '''
    df1 = pd.read_csv("country_k.csv")
    df2 = pd.read_csv("aT20200N10.csv")
    df1.drop(df1.tail(1).index, inplace=True)
    df2 = df2.loc[:, :"country_2"]
    df2 = df2.head(2000)
    result1 = pd.concat([df1["country_1"], df2["country_1"]], axis=1)
    result1.columns = ["Capital", "Productivity"]
    result1["Country"] = "country 1"
    Fig1 = sns.jointplot("Capital", "Productivity",data=result1, kind='hex')
    Fig1.fig.suptitle("An example of ergodic set with closed-formed solution", fontsize=16)
    Fig1.fig.subplots_adjust(top=0.95)
    return


def Figure2():
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

def Figure3(showcase_result):
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
    plt.ylabel("Mean Error($ln_{10}$)")
    ax = fig.add_subplot(1, 3, 3)
    plt.bar("Polynomial Degree", "Maximum Error", data=showcase_result, color=colors)
    plt.gca().invert_yaxis()
    plt.xlabel("Polynomial Degree")
    plt.ylabel("Maximum Error($ln_{10}$)")

    #build legend
    labels = ["Polynomial Degree "+str(i) for i in range(1, 6)]
    handles = [plt.Rectangle((0, 0), 1,1, color=colors[i]) for i in range(5)]
    plt.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), fancybox=True, shadow=True)
    plt.suptitle("Showcase results benchmarking", fontsize=16)
    plt.show()
    return

