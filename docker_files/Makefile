REPO=fischbachlab
NAME=nf-spades
TAG:=$(shell /bin/date +%Y%m%d%H%M%S)
REGISTRY=$(REPO)/$(NAME)

all: build push run

build:
	docker image build -t $(NAME):$(TAG) -t $(NAME):latest -f Dockerfile .
	docker image tag $(NAME):latest $(REGISTRY):latest
	docker image tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker image push $(REGISTRY):$(TAG)
	docker image push $(REGISTRY):latest
	echo "$(REPO)/$(NAME):$(TAG)" > LATEST

run:
	docker container run --rm $(NAME):latest
