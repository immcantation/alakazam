ARG REMOTE_USER=vscode
ARG REMOTE_UID=1000
ARG REMOTE_GID=1000

FROM immcantation/test:devel

LABEL maintainer="Susanna Marquez [susanna.marquez@yale.edu]" \
    description="Immcantation development container for testing"

ARG REMOTE_USER
ARG REMOTE_UID
ARG REMOTE_GID

# Install Docker
RUN dnf -y install dnf-plugins-core gh && \
    dnf-3 config-manager --add-repo https://download.docker.com/linux/fedora/docker-ce.repo && \
    dnf install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Install Act
RUN curl https://raw.githubusercontent.com/nektos/act/master/install.sh | bash -s -- -b /usr/local/bin

RUN getent group ${REMOTE_GID} || groupadd --non-unique --gid ${REMOTE_GID} ${REMOTE_USER}
RUN useradd --non-unique --uid ${REMOTE_UID} --gid ${REMOTE_GID} ${REMOTE_USER} \
    && echo "${REMOTE_USER} ALL=(root) NOPASSWD:ALL" > /etc/sudoers.d/${REMOTE_USER} \
    && chmod 0440 /etc/sudoers.d/${REMOTE_USER} \
    && echo "umask 000" >> /home/${REMOTE_USER}/.bashrc \
    && echo "umask 000" >> /root/.bashrc

# Add the user to the docker group
RUN usermod -aG docker ${REMOTE_USER}

USER ${REMOTE_USER}