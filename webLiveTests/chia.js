import {Selector} from 'testcafe';

fixture `Getting Started`
    .page `http://chia.systemsbiology.net`;

test('launch chiaa app', async function foo(t){

    const pageLoaded = Selector("#selectDatasetMenuLabel");
    await t
        .expect(pageLoaded.exists).ok()
    });
