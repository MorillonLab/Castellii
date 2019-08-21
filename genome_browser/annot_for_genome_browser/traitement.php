<?php 
  $twopoints = preg_match("#:#",$_POST['gene']);
  if($twopoints==0)
  {
    header('Location: images/'.strtoupper($_POST['gene']).'_'.$_POST['visu'].'.php');
  }
  else {
    if($twopoints==1) {
      $split = explode(":",$_POST['gene']);
      $decal = intval($_GET['decal']);
      $chr = strtolower($split[0]);
      $pos = max(1,ceil((floatval($split[1])/$decal)-1)*$decal);
      header('Location: images/'.$chr.'_'.$pos.'_'.$_POST['visu'].'.htm');
    }
    else{
      header('Location: 404.htm');
    }
  }
?>