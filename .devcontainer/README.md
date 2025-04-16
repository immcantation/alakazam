# Using .devcontainer for Development

This project includes **two `.devcontainer` configurations** to set up consistent development environments using [Visual Studio Code](https://code.visualstudio.com/) and [Docker](https://www.docker.com/):

1. **Testing Container**: Based on `immcantation/test:devel`. Smaller and with only the necessary tools for developing and testing this specific package.
2. **Development Container [Note: in preparation]**: Based on `immcantation/suite:devel`. Larger and with all Immcantation tools installed. Recommended if you are also working on other packages in the Immcantation suite or developing pipelines.


## Prerequisites

1. Install [Docker](https://docs.docker.com/get-docker/).
2. Install [Visual Studio Code](https://code.visualstudio.com/).
3. Install the [Dev Containers extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) for VS Code.
4. Set the environment variable GITHUB_TOKEN, with a [GitHub personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens). 

## Steps to Use the Dev Container

1. Clone this repository:
   ```bash
   git clone https://github.com/immcantation/alakazam.git
   cd alakazam
   ```

2. Open the repository in Visual Studio Code:
   ```bash
   code .
   ```

3. When prompted by VS Code, click **Reopen in Container** and select the container you want to use. This will build the development container and open the project inside it.

4. Start coding! The container includes all necessary dependencies for the project. To run tests, use the integrated terminal in VS Code. You can execute the `ci.yml` from the GitHub Local Actions plugin.

5. To stop the container, simply close VS Code or use the Docker extension to stop the running container.

## Additional Resources

- [Developing inside a Container (VS Code Docs)](https://code.visualstudio.com/docs/devcontainers/containers)
- [Docker Documentation](https://docs.docker.com/)
- [GitHub Local Actions](https://sanjulaganepola.github.io/github-local-actions-docs/)
