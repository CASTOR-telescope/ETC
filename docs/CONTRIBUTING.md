# Contributing

## Reporting Issues, Bugs, or Feature Requests

You can feel free to open a thread in [GitHub](https://github.com/CASTOR-telescope/ETC/discussions) for any discussions, or just directly [open an issue using one of our templates](https://github.com/CASTOR-telescope/ETC/issues/new/choose).

We ask that you create an issue for any features or changes you make so that Pull Requests can be linked accordingly. **Please create an issue prior to creating a pull request**. 

## Branch Strategy

The FORECASTOR team follows a [trunk-based development strategy](https://trunkbaseddevelopment.com/). This is enforced by repository settings on GitHub that will prevent direct pushes to the trunk branch.

![a picture of a scaled trunk-based development process, where each developer develops in their own branches](https://trunkbaseddevelopment.com/trunk1c.png)

**Please work in your own separate branch when contributing to FORECASTOR tools.** We recommend that you following our branch naming strategy of either `feat/<description>` or `issue/<issue-number>`

> To demonstrate the branch naming strategy with an example: Suppose I was working on a feature to add a "custom telescope" that implements issue #2 on GitHub. Two possible names I would use are `feat/custom-scope` or `issue/#2`