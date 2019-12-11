FROM rust:1.39.0-stretch as rust
RUN git clone https://github.com/hurwitzlab/trim-galore.git
RUN cd trim-galore/run_trim_galore && cargo build --release && cp target/release/run_trim_galore /usr/local/bin

FROM python:3.7.4-buster
RUN apt-get -y update
RUN apt-get install -y parallel wget git cutadapt default-jre
COPY --from=rust /usr/local/bin /usr/local/bin
COPY --from=rust /usr/local/lib /usr/local/lib
ENV LD_LIBRARY_PATH /usr/local/lib

# Install SRA Toolkit for fastq-dump
WORKDIR /app
RUN wget -q "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.8.zip" && unzip fastqc_v0.11.8.zip && ln -s /app/FastQC/fastqc /usr/local/bin/fastqc && chmod +x /usr/local/bin/fastqc && rm *.zip

RUN wget -q -O "trim_galore.tar.gz" "https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz" && tar -xvzf "trim_galore.tar.gz" && ln -s /app/TrimGalore-0.6.5/trim_galore /usr/local/bin/trim_galore && rm *.gz

# Run when the container launches
CMD ["run_trim_galore"]
