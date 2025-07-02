FROM rootproject/root:6.30.06-ubuntu22.04

# Install system dependencies for RapidSim
RUN apt-get update && apt-get install -y \
    git \
    cmake \
    make \
    libx11-dev \
    && rm -rf /var/lib/apt/lists/*

# Clone and build RapidSim
WORKDIR /opt
RUN git clone https://github.com/gcowan/RapidSim.git
WORKDIR /opt/RapidSim
RUN mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/opt/RapidSim && \
    make -j2 && \
    make install

# Set up RapidSim environment
ENV RAPIDSIM_ROOT=/opt/RapidSim
ENV PATH="${RAPIDSIM_ROOT}/build/src:$PATH"

# Copy project
WORKDIR /home/user
COPY . .

# Run ROOT scripts
WORKDIR /home/user/Setting
RUN root -l -b -q settingDaughters.cpp && \
    root -l -b -q settingMotherDst.cpp && \
    cp -r RapidSimConfigFiles/* /opt/RapidSim/config/smear/

# Make scripts executable
WORKDIR /home/user
RUN chmod +x run_pipeline.sh Simulation/simulation.sh

ENTRYPOINT ["./run_pipeline.sh"]