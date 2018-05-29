//------------------------------------------------------------------------------------------------------------------------
function setupMessageHandlers()
{
    Shiny.addCustomMessageHandler("roi",
                                  function(message) {
                                      console.log("roi message");
                                      console.log(message);
                                      window.igvBrowser.search(message.roi);
                                  }
                                 );

    Shiny.addCustomMessageHandler("displaySnps",

                                  function(message){
                                      console.log("displaySnps message");
                                      console.log(message);
                                      var trackName = "SNPs";
                                      var bedFileName = message.filename;
                                      var displayMode = "EXPANDED";
                                      var color = "red";
                                      var minValue = 0;
                                      var maxValue = 10;
                                      var trackHeight = 120;
                                      var url = window.location.href + bedFileName;
                                      var config = {format: "bedgraph",
                                                    name: trackName,
                                                    url: url,
                                                    min: minValue,
                                                    max: maxValue,
                                                    indexed: false,
                                                    displayMode: displayMode,
                                                    sourceType: "file",
                                                    color: color,
                                                    height: trackHeight,
                                                    type: "wig"};
                                      console.log("=== about to loadTrack");
                                      console.log(config)
                                      window.igvBrowser.loadTrack(config);
                                  })

    Shiny.addCustomMessageHandler("displayBedTrack", displayBedTrack);
    Shiny.addCustomMessageHandler("displayBedGraphTrack", displayBedGraphTrack);
    Shiny.addCustomMessageHandler("removeTrack", removeTrack);
    Shiny.addCustomMessageHandler("showGenomicRegion", showGenomicRegion);

} // setupMessageHandlers
//------------------------------------------------------------------------------------------------------------------------
function displayBedTrack(message)
{
  console.log("displayBedTrack, message:");
  console.log(message);

  var url = window.location.href + message.filename;

  var config = {format: "bed",
                name: message.trackName,
                url: url,
                indexed: false,
                displayMode: message.displayMode,
                sourceType: "file",
                color: message.color,
                height: message.trackHeight,
                type: "annotation"};
  console.log("=== about to loadTrack");
  console.log(config)
  window.igvBrowser.loadTrack(config);

} // displayBedTrack
//------------------------------------------------------------------------------------------------------------------------
function displayBedGraphTrack(message)
{
   console.log("displaySnps message");
   console.log(message);

   var trackName = message.trackName
   var filename = message.filename;
   var displayMode = "EXPANDED";
   var color = message.color;
   var minValue = message.minValue;
   var maxValue = message.maxValue;
   var trackHeight = message.trackHeight;

   var url = window.location.href + filename;

   var config = {format: "bedGraph",
                 name: trackName,
                 url: url,
                 min: 0,
                 max: 1,
                 indexed: false,
                 sourceType: "file",
                 color: color,
                 height: trackHeight,
                 type: "wig"};

   //console.log("=== about to loadTrack");
   //console.log(config)
   window.igvBrowser.loadTrack(config);

} // displayBedGraphTrack
//------------------------------------------------------------------------------------------------------------------------
function removeTrack(message)
{
  console.log("removeTrack, message:");
  console.log(message);

   var trackNames = message.trackName;

   if(typeof(trackNames) == "string")
      trackNames = [trackNames];

   var count = window.igvBrowser.trackViews.length;

   for(var i=(count-1); i >= 0; i--){
     var trackView = window.igvBrowser.trackViews[i];
     var trackViewName = trackView.track.name;
     var matched = trackNames.indexOf(trackViewName) >= 0;
     if (matched){
        window.igvBrowser.removeTrack(trackView.track);
        } // if matched
     } // for i

} // removeTrack
//------------------------------------------------------------------------------------------------------------------------
function showGenomicRegion(message)
{
  console.log("showGenomicRegion, message:");
  console.log(message);

  window.igvBrowser.search(message.region)

} // displayBedTrack
//------------------------------------------------------------------------------------------------------------------------
setTimeout(setupMessageHandlers, 5000);

