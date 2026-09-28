# How to contribute

The following has been adapted from Google's [open source code of conduct](https://opensource.google/documentation/reference/releasing/template/CODE_OF_CONDUCT). When in doubt, defer to that code of conduct



## Community Guidelines
### Be kind, empathetic, and helpful, especially to new contributors

### Resolve conflict constructively and peacefully

Healthy conflict, when resolved properly is vital to the success of projects like this where we want to make sure everyone's voice is heard, and credit for work is given appropriately. We especially encourage you to make your voice heard.
We believe however it needs to be resolved with high levels of respect and trust and *treating anyone with disrespect, aggression, or verbal abuse* is not okay. That being said, constructive conflict resolution requires high levels of effort
from all parties and we reserve the right to block/ban people (e.g if we're overwhelmed by low quality AI generated code/if someone acts in bad faith)

## Contribution guidelines

We welcome PRs and see the development section of the README to set up a development environment. 


### Volume

Contributions, especially from first time contributors should not be massive. A human will have to go through your code. We encourage you to open a draft PR and raise an issue if you want to work on a larger contribution, so that your priorities are in line
with the project and no work is wasted. In general if your first draft is AI generated, say so, explain how you've tried to verify behavior, what you expect it to do


### AI Generated code
This section is adapted from [Jax's guidelines](https://docs.jax.dev/en/latest/contributing.html). The underlying guideline for this section is that the maintainer should not have to do more work than what you have done yourself,
and you should not dump a massive amount of work on maintainers. We do expect to see AI generated code in commits and we do use tools like Claude Code

### Tests
In general any code that is contributed should 1. Not break any tests 2. Not introduce hidden regressions in end-use case that are not covered by tests. current AI is especially adept at convincing you that 2. will not happen, so be careful. 

The first check if you modify a test is ask why? and you should be able to explain that in your PR. Explain also what checks you've done to ensure that the code shows intended behavior

#### Responsiblity

You are responsible for every line of code you contribute and you should be able to explain what it does in context.

#### Communication
Do not use AI to speak for you. Any communication between you and a maintainer should solely be your words, barring edits. If the maintainers want to make a chatbot do work, we can do it ourselves, with more control over the outcome to boot.













