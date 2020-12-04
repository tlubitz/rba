function delete_session() {
    $.get('/simulator/clearsession', function(data) {
        console.log('Deleting session');
        window.location.href="http://127.0.0.1:8000/simulator";
    });
}

function simulate() {
    $.get('/simulator/simulate', function(data) {
        console.log('Simulating...');
        window.location.href="http://127.0.0.1:8000/simulator";
    });
}

function load_model(modelname) {
    var model_name_dict = {'modelname': modelname}

    // send the modelname to the view
    $.post('/simulator/loadmodel/', JSON.stringify(model_name_dict), function(data) {
        console.log('Choosing model from DB');
        window.location.href="http://127.0.0.1:8000/simulator";
    });

}