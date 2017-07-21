FROM ubuntu:16.04
RUN apt update
RUN apt install -y wget
RUN wget https://github.com/dzif/kd-tree-overlapper/releases/download/v1.0/kd -P /
RUN chmod a+x /kd

FROM bioboxes/biobox-minimal-base
RUN apt update && apt install -y patch
COPY --from=0 /kd  /usr/local/bin
COPY image/ /usr/local
RUN patch -p1 /usr/local/bin/validate_inputs.sh < /usr/local/validate.patch
ENV BIOBOX_EXEC execute_biobox.sh
ENV TASKFILE /usr/local/Taskfile
