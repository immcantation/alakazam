$(document).ready(function(){
    // Load this package announcements and insert into a new div
    $.get("announcements.html", function(data) {
        // Create a new div and add the fetched content
        const newDiv = $("<div>").html(data);
        
        // Insert the new div after #_1 and before #_2
        $("#_1").after(newDiv);
    }).fail(function() {
        console.error("Failed to load the message from the external file.");
    });    
});