ARG PIXI_VERSION=0.68.1

FROM ghcr.io/prefix-dev/pixi:${PIXI_VERSION} AS builder

WORKDIR /project
COPY pixi.toml pixi.lock run.py ./

RUN pixi install --locked
RUN pixi run build
RUN mv .pixi/envs/default /pixi-env

FROM debian:bookworm-slim

COPY --from=builder /pixi-env /pixi-env
ENV PATH="/pixi-env/bin:$PATH"

CMD ["/bin/bash"]
