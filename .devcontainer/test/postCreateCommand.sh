# #!/usr/bin/env bash

# # Install Act
# curl https://raw.githubusercontent.com/nektos/act/master/install.sh | sudo bash -s -- -b /usr/local/bin
# sudo dnf -y install dpkg

# # Install Docker
# sudo dnf -y install dnf-plugins-core gh
# sudo dnf-3 config-manager --add-repo https://download.docker.com/linux/fedora/docker-ce.repo
# sudo dnf install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

Rscript -e "dir.create(path = Sys.getenv('R_LIBS_USER'), showWarnings = FALSE, recursive = TRUE)"
Rscript -e 'install.packages(c("languageserver", "versions"))'

# Check if dependencies file exists and install
if [ -f "tests/setup/install_dep.R" ]; then
    Rscript tests/setup/install_dep.R
fi

echo -e "\n\033[1;32m============================================================"
echo -e " IMPORTANT: To install this R package, run in R console:"
echo -e "   devtools::install('.')"
echo -e "============================================================\033[0m\n"
