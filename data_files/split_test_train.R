split_test_train <- function(object, train_size){
  samp <- sample(nrow(object),floor(train_size*nrow(object)))
  train <- object[samp,]
  test <- object[-samp,]
  count <- 1

  #Repeat split 100x rejecting those with less than 5 presences in the test dataset and in the training dataset.
  while(sum(train[,1])< 5 | sum(test[,1])< 5){count <- count+1; if (count >100) break;  samp <- sample(nrow(object),floor(train_size*nrow(object)));  train <- object[samp,]; test <- object[-samp,]}
  
  return(list(samp,train,test,count))
  }