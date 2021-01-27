// determine dev/prod mode
if (window.location.href.substr(0,10) === 'http://127') { var mode = 'dev'; }
else { var mode = 'prod'; }

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
    $.get('/simulator/simulate', function(data) {
        console.log('Simulating...');
        if (mode == 'dev') { window.location.href="http://127.0.0.1:8000/simulator"; }
        else { window.location.href="https://timosan.pythonanywhere.com/simulator/"; }
    });
}
