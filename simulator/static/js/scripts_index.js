function delete_session() {
    $.get('/simulator/clearsession', function(data) {
        console.log('Deleting session');
        window.location.href="http://127.0.0.1:8000/simulator";
    });
}