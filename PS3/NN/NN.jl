using Flux
# using Flux: logitcrossentropy, normalise, onecold, onehotbatch
using Statistics: mean
using IterTools
using Parameters: @with_kw
using Tables


##

@with_kw mutable struct Args
    lr::Float64 = 0.5
    repeat::Int = 110
end
# args=Args()

##

include("ABKKExperiment.jl")

## %%
function get_processed_data(data,args)
    labels = vec(data."choice")
    features = Matrix(select(data, Not("choice")))'

    # Subract mean, divide by std dev for normed mean of 0 and std dev of 1.
    features = Flux.normalise(features, dims=2 )

    klasses = sort(unique(labels))
    onehot_labels = Flux.onehotbatch(labels, klasses)

    # Split into training and test sets, 2/3 for training, 1/3 for test.
    train_indices = [1:3:size(data,2) ; 2:3:size(data,2)]

    X_train = features[:,train_indices]
    y_train = onehot_labels[:,train_indices]

    X_test = features[:,3:3:size(data,1)]
    y_test = onehot_labels[:,3:3:size(data,1)]

    #repeat the data `args.repeat` times
    # train_data = Iterators.repeated([(X_train, y_train)],args.repeat)
    train_data = IterTools.ncycle([(X_train, y_train)],args.repeat)
    test_data = (X_test,y_test)

    return train_data, test_data, klasses
end


# train_data, test_data, klasses = get_processed_data(dnn,args)

## Accuracy Function
accuracy(x, y, model) = mean( Flux.onecold( model(x)) .== Flux.onecold(y) )

## Function to build confusion matrix
function confusion_matrix(X, y, model)
    args = Args()
    train_data, test_data, klasses = get_processed_data(dnn,args)
    ŷ = Flux.onehotbatch( Flux.onecold( model(X)) , 1:length(klasses))
    y * transpose(ŷ)
end
##
function train(model,args)
    # Initialize hyperparameter arguments
    # args = Args()

    #Loading processed data
    train_data, test_data, klasses = get_processed_data(dnn,args)

    # Declare model taking 4 features as inputs and outputting 3 probabiltiies,
    # one for each species of iris.
    # model = Chain(Dense(3177, 128, σ),
    #                 BatchNorm(128, relu),
    #                 Dense(128,64),
    #                 BatchNorm(64,relu),
    #                 Dense(64,32),
    #                 BatchNorm(32,relu),
    #                 Dense(32,10),
    #                 Dense(10,length(klasses)),
    #                 softmax)

    #Regularisation of Parameters
    # sqnorm(x) = sum(abs2, x)
    # Defining loss function to be used in training
    # For numerical stability, we use here logitcrossentropy
    loss(x, y) = Flux.logitcrossentropy(model(x), y)# + sum(sqnorm, Flux.params(model))

    # Training
    # Gradient descent optimiser with learning rate `args.lr`
    optimiser = Descent(args.lr)

    println("Starting training.")

    Flux.train!(loss, params(model), train_data, optimiser)

    x, y = test_data
    lss = loss(x,y)

    return test_data, lss

end

function test(model, test)
    # Testing model performance on test data
    X_test, y_test = test
    accuracy_score = accuracy(X_test, y_test, model)

    println("\nAccuracy: $accuracy_score")

    # Sanity check.
    # @assert accuracy_score > 0.8

    # To avoid confusion, here is the definition of a Confusion Matrix: https://en.wikipedia.org/wiki/Confusion_matrix
    println("\nConfusion Matrix:\n")
    display(confusion_matrix(X_test, y_test, model))

    cmat = confusion_matrix(X_test, y_test, model)

    return accuracy_score, cmat
end
## Calibrating
# cd(@__DIR__)

args = Args(lr=1.5,repeat=180)

# model = Chain(Dense(3177, 128, σ),
#                 BatchNorm(128, relu),
#                 Dense(128,64),
#                 BatchNorm(64,relu),
#                 Dense(64,32),
#                 BatchNorm(32,relu),
#                 Dense(32,10),
#                 Dense(10,length(klasses)),
#                 softmax)

model = Chain(Dense(3177,6))

## Testing

println("\n Args: learning rate = $(args.lr); Repetitions = $(args.repeat)")
println("\n Model: $model \n")
println("Starting training!")
test_data, lss = train(model,args)
acc, cmat = test(model, test_data)

## Saving results

CSV.write(tempdir1*"results.csv", DataFrame(Loss = lss, Accuracy = acc))
CSV.write(tempdir1*"cmat.csv", Tables.table(cmat))
