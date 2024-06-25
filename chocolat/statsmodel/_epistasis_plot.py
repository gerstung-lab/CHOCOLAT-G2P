import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np 
import statsmodels.api as sm
import warnings
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning


# Define the model
class BinomialRegressionModel(nn.Module):
    """
    A binomial regression model implemented in torch.

    Parameters
    ----------
    num_samples : int
        Number of samples.
    input_features : int
        Number of input features.
    output_features : int
        Number of output features.
    """

    def __init__(self, num_samples, input_features, output_features):
        super(BinomialRegressionModel, self).__init__()
        self.weights = nn.Parameter(torch.randn(num_samples, input_features, output_features))
    
    def forward(self, x):
        """
        Forward pass for the binomial regression model.

        Parameters
        ----------
        x : torch.Tensor
            Input tensor.

        Returns
        -------
        torch.Tensor
            Output tensor after applying the model weights.
        """

        weights = torch.bmm(x, self.weights)
        return weights

        
def plot_epistatic_effect_torch(sim_genotypes, a,b, plasmids_ordered_list, ax=None):
    """
    Plots the epistatic comparison of simulated genotypes between real and expected data, using torch model.

    Parameters
    ----------
    sim_genotypes : numpy.ndarray
        Simulated genotypes (n_plasmids, n_samples[at least 2000], n_nodules).
    a : int
        Index of the first plasmid.
    b : int
        Index of the second plasmid.
    plasmids_ordered_list : list of str
        List of ordered plasmids (names).
    ax : matplotlib.axes.Axes, optional
        Axis object to draw the plot onto, otherwise creates new figure.

    Returns
    -------
    None
    """

    dat_real = np.array([np.apply_along_axis(lambda x: list(x) == [0,0], 0, sim_genotypes[[a,b],:2000]),
                          np.apply_along_axis(lambda x: list(x) == [1,0], 0, sim_genotypes[[a,b],:2000]),
                          np.apply_along_axis(lambda x: list(x) == [0,1], 0, sim_genotypes[[a,b],:2000]),
                          np.apply_along_axis(lambda x: list(x) == [1,1], 0, sim_genotypes[[a,b],:2000])]).mean(-1).T

    # Define the dataset
    X = torch.tile(torch.tensor([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 1]], dtype=torch.float), (2000,1,1))
    X_al = torch.tile(torch.tensor([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], dtype=torch.float), (2000,1,1))

    y = torch.tensor(dat_real, dtype=torch.float) 


    model = BinomialRegressionModel(2000,3,1)

    # Loss and Optimizer
    criterion = nn.CrossEntropyLoss()  # Mean Squared Error Loss
    optimizer = optim.Adam(model.parameters(), lr=0.01)

    # Training the model
    num_epochs = 5000
    for epoch in range(num_epochs):
        # Forward pass: Compute predicted y by passing x to the model
        y_pred = model(X)

        # Compute and print loss
        loss = criterion(y_pred.squeeze(), y)

        # Zero gradients, perform a backward pass, and update the weights.
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        
    dat_real = torch.softmax(model(X), dim=1)[:,:,0].detach().numpy()
    dat_expected = torch.softmax(model(X_al), dim=1)[:,:,0].detach().numpy()

    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(4, 3))
        
    for i in range(4):
        ax.bar(np.arange(4)[i]-0.5, np.median(dat_real, 0)[i], width=0.2, color='tomato')
        ax.bar(np.arange(4)[i]-0.25, np.median(dat_expected, 0)[i], width=0.2, color='dodgerblue', alpha=0.3)

        ax.plot([np.arange(4)[i]-0.5, np.arange(4)[i]-0.5,], [np.percentile(dat_real, 2.5, 0)[i], np.percentile(dat_real, 97.5, 0)[i]], c='k', lw=1);
        ax.plot([np.arange(4)[i]-0.25, np.arange(4)[i]-0.25,], [np.percentile(dat_expected, 2.5, 0)[i], np.percentile(dat_expected, 97.5, 0)[i]], c='k', lw=1);

    ax.set_xticks(np.arange(4)-0.375)
    ax.set_xticklabels([f'{plasmids_ordered_list[::-1][a].split("-")[0]}:NO;{plasmids_ordered_list[::-1][b].split("-")[0]}:NO',
                              f'{plasmids_ordered_list[::-1][a].split("-")[0]}:YES;{plasmids_ordered_list[::-1][b].split("-")[0]}:NO',
                              f'{plasmids_ordered_list[::-1][a].split("-")[0]}:NO;{plasmids_ordered_list[::-1][b].split("-")[0]}:YES',
                              f'{plasmids_ordered_list[::-1][a].split("-")[0]}:YES;{plasmids_ordered_list[::-1][b].split("-")[0]}:YES'], rotation=90)



    ax.set_ylim(0,0.7);
    ax.spines[['right', 'top']].set_visible(False)
    

def plot_epistatic_effect_glm(sim_genotypes, a, b, plasmids_ordered_list, ax=None, n=5000):
    """
    Plots the epistatic effects of genotypes using a Generalized Linear Model (GLM).

    Parameters
    ----------
    sim_genotypes : numpy.ndarray
        Simulated genotypes (n_plasmids, n_samples, n_nodules).
    a : int
        Index of the first plasmid.
    b : int
        Index of the second plasmid.
    plasmids_ordered_list : list of str
        List of ordered plasmids.
    ax : matplotlib.axes.Axes, optional
        Axis object to draw the plot onto, otherwise uses the current axis.
    n : int, optional
        Number of samples to use. Default is 5000.

    Returns
    -------
    None
    """

    y = np.array([np.apply_along_axis(lambda x: list(x) == [0,0], 0, sim_genotypes[[a,b],:n]),
          np.apply_along_axis(lambda x: list(x) == [1,0], 0, sim_genotypes[[a,b],:n]),
          np.apply_along_axis(lambda x: list(x) == [0,1], 0, sim_genotypes[[a,b],:n]),
          np.apply_along_axis(lambda x: list(x) == [1,1], 0, sim_genotypes[[a,b],:n])]).mean(-1)

    ysh = y[:,0]

    x = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])

    # Create an interaction term
    x_interaction = x[:, 0] * x[:, 1]
    x_interaction = x_interaction.reshape(-1, 1)

    # Add the interaction term to your features
    x = np.concatenate([x, x_interaction], axis=1)
    x = sm.add_constant(x)
    
    warnings.filterwarnings("ignore", category=PerfectSeparationWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in scalar divide")
    warnings.filterwarnings("ignore", category=RuntimeWarning, message="invalid value encountered in scalar divide")

    params = []
    for i in range(n):
        model = sm.GLM(y[:,i], x, family=sm.families.Binomial())
        result = model.fit()
        params.append(result.params)
        
    params = np.array(params).T
    
    rel_space = np.array([1 / (1 + np.exp(-(params[0]))),
                        1 / (1 + np.exp(-(params[0] + params[1]))),
                        1 / (1 + np.exp(-(params[0] + params[2]))),
                        1 / (1 + np.exp(-(params[0] + params[1] + params[2] + params[3]))),
                        1 / (1 + np.exp(-(params[0] + params[1] + params[2])))])
    
    from matplotlib.markers import MarkerStyle

    dat_real = rel_space[[0,2,1,4]].T
    dat_expected = y.T

    if ax is None:
        fig, ax = plt.subplots(1,1,figsize=(4, 3))
    for i in range(4):
        ax.bar(np.arange(4)[i]-0.5, np.median(dat_real, 0)[i], width=0.2, color='tomato')
        ax.bar(np.arange(4)[i]-0.25, np.median(dat_expected, 0)[i], width=0.2, color='dodgerblue', alpha=0.3)

        ax.plot([np.arange(4)[i]-0.5, np.arange(4)[i]-0.5,], [np.percentile(dat_real, 2.5, 0)[i], np.percentile(dat_real, 97.5, 0)[i]], c='k', lw=1);
        ax.plot([np.arange(4)[i]-0.25, np.arange(4)[i]-0.25,], [np.percentile(dat_expected, 2.5, 0)[i], np.percentile(dat_expected, 97.5, 0)[i]], c='k', lw=1);


    ax.scatter(np.arange(4)[-1], np.median(y.T,0)[-1], c='w', s=15, zorder=3, edgecolor='k')

    ax.scatter(np.arange(4)[-1], np.median(rel_space[[0,2,1,4]].T,0)[-1], c='white', s=15, zorder=3, edgecolor='k')

    pv = np.minimum(np.mean(dat_real[:,-1] >= dat_expected[:,-1]), np.mean(dat_real[:,-1] < dat_expected[:,-1]))
    ax.set_title(f'{plasmids_ordered_list[::-1][a]}/{plasmids_ordered_list[::-1][b]} ; pval = {pv}')

    ax.arrow(3, np.median(dat_expected,0)[-1], 0,
              np.median(dat_real,0)[-1]-np.median(dat_expected,0)[-1],
              head_width=0.1, head_length=0.03, color='k', length_includes_head=True)



    ax.set_xticks(np.arange(4)-0.375)
    ax.set_xticklabels([f'{plasmids_ordered_list[::-1][a].split("-")[0]}:NO;{plasmids_ordered_list[::-1][b].split("-")[0]}:NO',
                              f'{plasmids_ordered_list[::-1][a].split("-")[0]}:NO;{plasmids_ordered_list[::-1][b].split("-")[0]}:YES',
                              f'{plasmids_ordered_list[::-1][a].split("-")[0]}:YES;{plasmids_ordered_list[::-1][b].split("-")[0]}:NO',
                              f'{plasmids_ordered_list[::-1][a].split("-")[0]}:YES;{plasmids_ordered_list[::-1][b].split("-")[0]}:YES'], rotation=90)



    ax.set_ylim(0,0.7);
    ax.spines[['right', 'top']].set_visible(False)

