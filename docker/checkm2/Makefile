IMAGE_NAME = checkm2
VERSION = 1.1.0

REPOS = \
	ghcr.io/broadinstitute/hvp-lr \
	us-central1-docker.pkg.dev/broad-hvp-dasc/hvp-longread-containers

################################################################################
################################################################################
#       ____        _ _             ____  _       _
#      | __ )  ___ (_) | ___ _ __  |  _ \| | __ _| |_ ___
#      |  _ \ / _ \| | |/ _ \ '__| | |_) | |/ _` | __/ _ \
#      | |_) | (_) | | |  __/ |    |  __/| | (_| | ||  __/
#      |____/ \___/|_|_|\___|_|    |_|   |_|\__,_|\__\___|
#
################################################################################
################################################################################

########################################
# Buildx settings
########################################

# Multi-arch platforms for `push`. `build` stays single-arch for fast local dev.
PLATFORMS = linux/amd64,linux/arm64

# Force Docker v2.2 manifest list (application/vnd.docker.distribution.manifest.list.v2+json).
# Some downstream consumers reject OCI mediatypes.
# provenance + sbom MUST be off — both attestations force OCI mediatypes.
BUILDX_ATTEST_OFF = --provenance=false --sbom=false
MANIFEST_OUTPUT   = type=image,oci-mediatypes=false,push=true

########################################
# Set up image tags
########################################

TAGS = $(VERSION) latest main
TARGETS  = $(foreach repo,$(REPOS),$(foreach tag,$(TAGS),$(repo)/$(IMAGE_NAME):$(tag)))
TAG_ARGS = $(foreach t,$(TARGETS),-t $(t))

########################################
# Build targets
########################################

# Reminder: ORDER MATTERS HERE - PUT `all` first!
all: build push

# Stage the repo-root LICENSE into the build context so Dockerfile's
# `COPY LICENSE` works. Regenerated whenever the root LICENSE changes.
# Ignored by git (see root .gitignore).
LICENSE: ../../LICENSE
	cp ../../LICENSE LICENSE

# Local single-arch build. Loads image into the local docker daemon as
# $(IMAGE_NAME):$(VERSION) and :latest. Fast — no cross-arch QEMU.
build: LICENSE
	docker buildx build --load --platform linux/amd64 \
	    $(BUILDX_ATTEST_OFF) \
	    -t $(IMAGE_NAME):$(VERSION) -t $(IMAGE_NAME):latest .

build_no_cache: LICENSE
	docker buildx build --no-cache --load --platform linux/amd64 \
	    $(BUILDX_ATTEST_OFF) \
	    -t $(IMAGE_NAME):$(VERSION) -t $(IMAGE_NAME):latest .

# Internal name for the multi-arch buildx builder used by `push`.
BUILDER = hvp-builder

# Multi-arch build + push to every entry in REPOS as Docker v2.2 manifest list.
# Requires `docker login` to each registry in REPOS. The multi-arch builder is
# created on demand if it doesn't already exist.
# Does NOT depend on `build` — buildx builds and pushes in one shot from source.
push: LICENSE
	@docker buildx inspect $(BUILDER) >/dev/null 2>&1 || \
	    docker buildx create --name $(BUILDER) --driver docker-container --bootstrap >/dev/null
	docker buildx build --builder $(BUILDER) --platform $(PLATFORMS) \
	    $(BUILDX_ATTEST_OFF) \
	    --output $(MANIFEST_OUTPUT) \
	    $(TAG_ARGS) .

clean:
	rm -f LICENSE

.PHONY: all build build_no_cache push clean
