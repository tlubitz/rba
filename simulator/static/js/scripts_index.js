// determine dev/prod mode
if (window.location.href.substr(0,10) === 'http://127') { var mode = 'dev'; }
else { var mode = 'prod'; }

var coll = document.getElementsByClassName("collapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    if (content.style.maxHeight){
      content.style.maxHeight = null;
    } else {
      content.style.maxHeight = content.scrollHeight + "px";
    }
  });
}


function delete_session() {
    $.get('/simulator/clearsession', function(data) {
        console.log('Deleting session');
        if (mode == 'dev') { window.location.href="http://127.0.0.1:8000/simulator"; }
        else { window.location.href="https://timosan.pythonanywhere.com/simulator/"; }
    });
}


function load_model(modelname) {
    var model_name_dict = {'modelname': modelname}

    // send the modelname to the view
    $.post('/simulator/loadmodel/', JSON.stringify(model_name_dict), function(data) {
        console.log('Choosing model from DB');
        if (mode == 'dev') { window.location.href="http://127.0.0.1:8000/simulator"; }
        else { window.location.href="https://timosan.pythonanywhere.com/simulator/"; }
    });
}


function simulate() {
    var parameters_species = {}
    var parameters = {}
    var species = {}

    // first collect available parameters
    var parameter_container = document.querySelector('.parameters');
    var p_inputs = parameter_container.querySelectorAll('input');
    for (var i=0; i<p_inputs.length; i++) {
        var v = p_inputs[i].value;
        if (v == '') { v = p_inputs[i].placeholder; }
        parameters[p_inputs[i].id] = v;
    }
    parameters_species['parameters'] = parameters;

    // then the species concentrations
    var species_container = document.querySelector('.specieses');
    var s_inputs = species_container.querySelectorAll('input');
    for (var i=0; i<s_inputs.length; i++) {
        var v = s_inputs[i].value;
        if (v == '') { v = s_inputs[i].placeholder; }
        species[s_inputs[i].id] = v;
    }
    parameters_species['species'] = species;

    $.post('/simulator/simulate/', JSON.stringify(parameters_species), function(data) {
        console.log('Simulating...');
        if (mode == 'dev') { window.location.href="http://127.0.0.1:8000/simulator"; }
        else { window.location.href="https://timosan.pythonanywhere.com/simulator/"; }
    });
}


function undo_last() {
  $.get('/simulator/undolast/', function(data) {
      console.log('Undo last step');
      if (mode == 'dev') { window.location.href="http://127.0.0.1:8000/simulator"; }
      else { window.location.href="https://timosan.pythonanywhere.com/simulator/"; }
  });
}


function plot_parameter() {
  var select_box = document.querySelector('#select_parameter');
  parameter = {};
  parameter['plot_parameter'] = select_box.value;

  $.post('/simulator/plot/', JSON.stringify(parameter), function(data) {
    console.log('Plotting...');
    if (mode == 'dev') { window.location.href="http://127.0.0.1:8000/simulator#plot-path"; }
    else { window.location.href="https://timosan.pythonanywhere.com/simulator#plot-path"; }
  });
}