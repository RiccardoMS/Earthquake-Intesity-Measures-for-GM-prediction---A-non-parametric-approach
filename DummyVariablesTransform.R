DummyVariablesTransform<-function (lungh)
{ for (i in 1:lungh){
  if(EC8Group=="A"){
    DummyA=1
    break()}
  else
    DUmmyA=0
  
  if(EC8Group=="A*"){
    DummyAStar=1
    break()}
  else
      DummyAStar=0
  
  if(EC8Group=="B"){
    DummyB=1
    break()}
  else
    DummyB=0
  
  if(EC8Group=="B*"){
    DummyBStar=1
    break()}
  else
    DummyBStar=0
  
  if(EC8Group=="C"){
    DummyC=1
    break()}
  else
    DummyC=0
  
  if(EC8Group=="C*"){
    DummyCStar=1
    break()}
  else
    DummyCStar=0
  
  if(EC8Group=="E"){
    DummyE=1
    break()}
  else
    DummyE=0
}
DATA<-cbind(DATA,DummyA,DummyAStar,DummyB,DummyBStar,DummyC,DummyCStar,DummyE) 
}