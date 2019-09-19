# Package Guide

## Installation

GenomeGraphs is made available to install through BioJulia's package registry.

Julia by default only watches the "General" package registry, so before you start, you should add the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the prompt change from `julia>` to `pkg>`), then enter the following command:

```
pkg> registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install FASTX from the julia REPL. Press ] to enter pkg mode again, and enter the following:

```
pkg> add GenomeGraphs
```

If you are interested in the cutting edge of the development, please check out the master branch to try new features before release.

## Creating a genome graph WorkSpace

You can create an empty genome graph workspace with the empty constructor:

```@example
using GenomeGraphs

ws = WorkSpace()
```

Obviously this workspace is quite useless on its own!

You need a genome graph and information to project onto it before you can
do any exploration or analysis. There are a few ways to get your first graph.

### Constructing a de-bruijn graph

