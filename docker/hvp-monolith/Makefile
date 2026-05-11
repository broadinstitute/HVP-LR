IMAGE_NAME = hvp-monolith
VERSION = 0.0.1

REPOS = \
	us.gcr.io/broadinstitute \
	ghcr.io/broadinstitute

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
# Set up image names and build function:
########################################

TAGS = $(VERSION) latest main
TARGETS = $(foreach repo, $(REPOS), $(foreach tag, $(TAGS), $(repo)/$(IMAGE_NAME):$(tag)))

define tag_and_push
	docker tag $(IMAGE_NAME) $(1)
	docker push $(1)

endef

########################################
# Define build targets:
########################################

# Reminder: ORDER MATTERS HERE - PUT `all` first!
all: build push

build:
	docker build -t $(IMAGE_NAME):$(version) .

build_no_cache:
	docker build --no-cache -t $(IMAGE_NAME):$(version) .

push: build
	$(foreach target, $(TARGETS), $(call tag_and_push, $(target)))

.PHONY: all build push

